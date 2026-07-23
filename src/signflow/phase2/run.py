"""Phase 2 orchestration for expression-only processing and DE."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from signflow.config import load_config
from signflow.provenance import DirtyRepositoryError, git_state, utc_now_path_component, write_run_manifest
from signflow.resources import PACKAGE_VERSION_NAMES

from . import PHASE2_META_ROOT, PHASE2_OUTPUT_ROOT, Phase2Result, StopAndReportError, ensure_directory
from .acquire import acquire_geo_series
from .activity import infer_activities
from .bootstrap_signature import bootstrap_meta_signatures
from .concordance import summarize_concordance
from .de_limma import limma_versions, run_limma_differential_expression, write_de_plots
from .meta_analysis import meta_analyze_signatures
from .phenotype import (
    classify_phenotypes,
    phenotype_mapping_summary,
    select_gse11783_primary_samples,
    select_uc_only,
)
from .process import process_expression_study
from .qc import write_qc_outputs

TARGET_PRIMARY = {
    "GSE11783": {"platform": "GPL570", "expected_case": 5, "expected_control": 6},
    "GSE57560": {"platform": "GPL16699", "expected_case": 13, "expected_control": 3},
}
TARGET_CONTROL = {"GSE238208": {}}
BENCHMARK = {"GSE75214": {"platform": "GPL6244"}, "GSE87466": {"platform": "GPL13158"}}
COMBAT_HEADER = (
    "SENSITIVITY ONLY — pooled ComBat-like study mean-centering. "
    "Study×condition confounding: GSE11783 6:5 vs GSE57560 3:13. Not primary."
)


def _count_groups(phenotypes: pd.DataFrame) -> dict[str, int]:
    return {str(key): int(value) for key, value in phenotypes["group"].value_counts().to_dict().items()}


def _assert_expected_counts(accession: str, phenotypes: pd.DataFrame, expected_case: int, expected_control: int) -> None:
    counts = _count_groups(phenotypes)
    if counts.get("case", 0) != expected_case or counts.get("control", 0) != expected_control:
        raise StopAndReportError(
            f"{accession} phenotype counts drifted from DATA_AVAILABILITY.md / frozen plan",
            accession=accession,
            reason="count-drift",
            details={"observed": counts, "expected": {"case": expected_case, "control": expected_control}},
        )


def _study_output_dir(accession: str) -> Path:
    return ensure_directory(PHASE2_OUTPUT_ROOT / accession)


def _save_table(frame: pd.DataFrame, path: str | Path, *, index: bool = False, header_comment: str | None = None) -> str:
    path = Path(path)
    ensure_directory(path.parent)
    if header_comment:
        path.write_text("# " + header_comment + "\n")
        frame.to_csv(path, sep="\t", index=index, mode="a")
    else:
        frame.to_csv(path, sep="\t", index=index)
    return str(path)


def _gse11783_variant_phenotypes(phenotypes_all: pd.DataFrame, mode: str) -> pd.DataFrame:
    if mode == "non-ulcer":
        selected, _ = select_gse11783_primary_samples(phenotypes_all, primary_site="non-ulcer")
        return selected
    if mode == "ulcer":
        selected, _ = select_gse11783_primary_samples(phenotypes_all, primary_site="ulcer")
        return selected
    if mode == "patient-mean":
        return phenotypes_all[phenotypes_all["group"].isin(["case", "control"])].reset_index(drop=True)
    raise ValueError(mode)


def _patient_mean_expression(expression: pd.DataFrame, phenotypes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    controls = phenotypes[phenotypes["group"] == "control"]
    cases = phenotypes[phenotypes["group"] == "case"]
    columns = {}
    rows = []
    for sample_id in controls["sample_id"]:
        columns[sample_id] = expression[sample_id]
        rows.append({"sample_id": sample_id, "group": "control", "patient_id": None, "site": "healthy"})
    for patient_id, group in cases.groupby("patient_id"):
        sample_ids = group["sample_id"].tolist()
        col = f"patient_{patient_id}"
        columns[col] = expression[sample_ids].mean(axis=1)
        rows.append({"sample_id": col, "group": "case", "patient_id": patient_id, "site": "patient-mean"})
    return pd.DataFrame(columns), pd.DataFrame(rows)


def _rank_combine(study_results: list[pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for result in study_results:
        frame = result[["gene_symbol", "effect"]].copy()
        frame["rank"] = frame["effect"].rank(method="average")
        study = str(result["study"].iloc[0])
        frames.append(frame[["gene_symbol", "rank"]].rename(columns={"rank": study}))
    merged = frames[0]
    for frame in frames[1:]:
        merged = merged.merge(frame, on="gene_symbol", how="inner")
    study_cols = [c for c in merged.columns if c != "gene_symbol"]
    merged["mean_rank"] = merged[study_cols].mean(axis=1)
    return merged.sort_values("mean_rank", ascending=False).reset_index(drop=True)


def _combat_like_pooled(
    study_expression: dict[str, pd.DataFrame],
    study_phenotypes: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    from .de_limma import standardized_effects

    centered_parts = []
    pheno_parts = []
    for study, expr in study_expression.items():
        centered = expr.sub(expr.mean(axis=1), axis=0)
        centered_parts.append(centered)
        pheno = study_phenotypes[study][["sample_id", "group"]].copy()
        pheno["study"] = study
        pheno_parts.append(pheno)
    genes = sorted(set.intersection(*(set(part.index) for part in centered_parts)))
    pooled = pd.concat([part.loc[genes] for part in centered_parts], axis=1)
    pheno = pd.concat(pheno_parts, ignore_index=True)
    effects = standardized_effects(pooled, pheno)
    effects["analysis"] = "combat_like_sensitivity"
    return effects


def _process_target_primary(
    accession: str,
    acquisition: dict[str, object],
    result: Phase2Result,
    *,
    gse11783_mode: str = "non-ulcer",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    phenotypes_all = classify_phenotypes(accession, acquisition["series_matrix_path"])
    study_dir = _study_output_dir(accession)
    result.notes.append(f"{accession} [{gse11783_mode}]: {phenotype_mapping_summary(phenotypes_all)}")

    if accession == "GSE11783":
        phenotypes = _gse11783_variant_phenotypes(phenotypes_all, gse11783_mode)
        if gse11783_mode == "non-ulcer":
            _, patient_map = select_gse11783_primary_samples(phenotypes_all, primary_site="non-ulcer")
            result.files_created.append(_save_table(patient_map, study_dir / "patient_sample_map.tsv"))
            _assert_expected_counts(
                accession,
                phenotypes,
                TARGET_PRIMARY[accession]["expected_case"],
                TARGET_PRIMARY[accession]["expected_control"],
            )
            result.per_accession_counts[accession] = _count_groups(phenotypes)
        else:
            result.notes.append(f"GSE11783 sensitivity mode={gse11783_mode}; primary count assert skipped.")
    else:
        phenotypes = phenotypes_all[phenotypes_all["group"].isin(["case", "control"])].reset_index(drop=True)
        _assert_expected_counts(
            accession,
            phenotypes,
            TARGET_PRIMARY[accession]["expected_case"],
            TARGET_PRIMARY[accession]["expected_control"],
        )
        result.per_accession_counts[accession] = _count_groups(phenotypes)

    pheno_name = "phenotype_table.tsv" if gse11783_mode == "non-ulcer" or accession != "GSE11783" else f"phenotype_table_{gse11783_mode}.tsv"
    result.files_created.append(_save_table(phenotypes, study_dir / pheno_name))

    sample_ids = phenotypes["sample_id"].tolist()
    expression, mapping_qc, processing_meta = process_expression_study(
        accession,
        TARGET_PRIMARY[accession]["platform"],
        acquisition["series_matrix_path"],
        sample_ids,
    )
    result.notes.append(f"{accession} processing: {processing_meta}")
    qc_name = "mapping_qc.tsv" if gse11783_mode == "non-ulcer" or accession != "GSE11783" else f"mapping_qc_{gse11783_mode}.tsv"
    result.files_created.append(_save_table(mapping_qc, study_dir / qc_name))

    if accession == "GSE11783" and gse11783_mode == "patient-mean":
        expression, phenotypes = _patient_mean_expression(expression, phenotypes)

    expression_path = study_dir / (
        "expression_gene_level.tsv" if gse11783_mode == "non-ulcer" or accession != "GSE11783" else f"expression_gene_level_{gse11783_mode}.tsv"
    )
    expression.to_csv(expression_path, sep="\t")
    result.files_created.append(str(expression_path))

    de_table = run_limma_differential_expression(expression, phenotypes, accession=accession)
    de_table["study"] = accession
    de_name = "de_limma.tsv" if gse11783_mode == "non-ulcer" or accession != "GSE11783" else f"de_limma_{gse11783_mode}.tsv"
    de_path = study_dir / de_name
    de_table.to_csv(de_path, sep="\t", index=False)
    result.files_created.append(str(de_path))
    result.files_created.extend(write_de_plots(de_table, study_dir / "plots", accession=f"{accession}_{gse11783_mode}"))

    if gse11783_mode == "non-ulcer" or accession != "GSE11783":
        result.files_created.extend(write_qc_outputs(expression, phenotypes, accession, study_dir / "qc"))
        try:
            result.files_created.extend(list(infer_activities(expression, study_dir / "activity").values()))
        except StopAndReportError as exc:
            result.stop_conditions.append(exc.to_record())
            result.notes.append(f"{accession}: activity inference stopped — {exc}")
    return de_table, expression, phenotypes


def _process_target_control(accession: str, acquisition: dict[str, object], result: Phase2Result) -> None:
    phenotypes = classify_phenotypes(accession, acquisition["series_matrix_path"])
    study_dir = _study_output_dir(accession)
    result.per_accession_counts[accession] = _count_groups(phenotypes)
    result.files_created.append(_save_table(phenotypes, study_dir / "phenotype_table.tsv"))
    result.notes.append(
        f"{accession}: acquired and phenotyped only. Inflammation-specificity control; "
        "forbidden as IC-vs-healthy DE."
    )


def _process_benchmark(accession: str, acquisition: dict[str, object], result: Phase2Result) -> None:
    phenotypes = classify_phenotypes(accession, acquisition["series_matrix_path"])
    if accession == "GSE75214":
        phenotypes, excluded = select_uc_only(phenotypes)
        result.notes.append(
            f"{accession}: excluded {excluded['excluded']} non-primary samples "
            "(CD/Crohn's/inactive UC/ileum). Frozen contrast: UC active colon vs control colon."
        )
    study_dir = _study_output_dir(accession)
    result.per_accession_counts[accession] = _count_groups(phenotypes)
    result.files_created.append(_save_table(phenotypes, study_dir / "phenotype_table.tsv"))


def run_phase2(*, track: str, allow_dirty: bool = False, n_bootstrap: int = 1000) -> Phase2Result:
    """Run Phase 2 for TARGET or BENCHMARK cohorts."""

    track = track.lower()
    if track not in {"target", "benchmark"}:
        raise ValueError("track must be 'target' or 'benchmark'")

    git = git_state(".")
    if git["dirty"] and not allow_dirty:
        raise DirtyRepositoryError("Repository is dirty; re-run with --allow-dirty to continue.")

    result = Phase2Result(track=track)
    downloads: list[dict[str, object]] = []
    config = load_config("configs/default.yaml")
    try:
        if track == "target":
            result.notes.append(f"limma/R versions: {limma_versions()}")
            study_results = []
            study_expression = {}
            study_phenotypes = {}
            acquisitions = {}
            for accession in TARGET_PRIMARY:
                acquisitions[accession] = acquire_geo_series(accession)
                downloads.extend(acquisitions[accession]["download_records"])
                de_table, expression, phenotypes = _process_target_primary(
                    accession, acquisitions[accession], result, gse11783_mode="non-ulcer"
                )
                study_results.append(de_table)
                study_expression[accession] = expression
                study_phenotypes[accession] = phenotypes

            for mode in ("ulcer", "patient-mean"):
                de_sens, _, _ = _process_target_primary(
                    "GSE11783", acquisitions["GSE11783"], result, gse11783_mode=mode
                )
                sens_meta = meta_analyze_signatures([de_sens, study_results[1]])
                sens_path = ensure_directory(PHASE2_META_ROOT / "sensitivity") / f"signature_meta_gse11783_{mode}.tsv"
                result.files_created.append(
                    _save_table(sens_meta, sens_path, header_comment=f"SENSITIVITY ONLY — GSE11783 {mode}")
                )

            meta = meta_analyze_signatures(study_results)
            meta_path = ensure_directory(PHASE2_META_ROOT) / "signature_meta.tsv"
            meta.to_csv(meta_path, sep="\t", index=False)
            result.files_created.append(str(meta_path))
            lost_a = len(set(study_results[0]["gene_symbol"]) - set(meta["gene_symbol"]))
            lost_b = len(set(study_results[1]["gene_symbol"]) - set(meta["gene_symbol"]))
            result.notes.append(
                f"Primary meta-analysis: DerSimonian-Laird on Hedges' g; "
                f"intersection genes={len(meta)}; lost GSE11783-only={lost_a}; lost GSE57560-only={lost_b}"
            )
            result.concordance = summarize_concordance(study_results[0], study_results[1])
            result.concordance["i2_distribution"] = {
                "median": float(meta["i2"].median()),
                "mean": float(meta["i2"].mean()),
                "q75": float(meta["i2"].quantile(0.75)),
            }
            result.notes.append(str(result.concordance.get("plain_statement", "")))

            rank_path = ensure_directory(PHASE2_META_ROOT / "sensitivity") / "rank_combine.tsv"
            result.files_created.append(
                _save_table(_rank_combine(study_results), rank_path, header_comment="SENSITIVITY ONLY — rank-based combination")
            )
            combat = _combat_like_pooled(study_expression, study_phenotypes)
            combat_path = ensure_directory(PHASE2_META_ROOT / "sensitivity") / "combat_like_effects.tsv"
            result.files_created.append(_save_table(combat, combat_path, header_comment=COMBAT_HEADER))

            meta_expr = pd.DataFrame({"meta": meta.set_index("gene_symbol")["pooled_effect"]})
            try:
                result.files_created.extend(list(infer_activities(meta_expr, PHASE2_META_ROOT / "activity_meta").values()))
            except StopAndReportError as exc:
                result.stop_conditions.append(exc.to_record())

            boot = bootstrap_meta_signatures(
                study_expression,
                study_phenotypes,
                n_bootstrap=n_bootstrap,
                seed=0,
                output_dir=PHASE2_META_ROOT / "bootstrap",
            )
            result.files_created.extend([boot["summary_path"], boot["index_path"]])
            result.notes.append(
                f"Bootstrap: n={boot['n_bootstrap']}, stratified preserving group sizes; "
                f"degenerate(sum over studies)={boot['degenerate_resamples_sum_over_studies']}; "
                "propagates Hedges'g→DL meta (not limma each resample)."
            )

            for accession in TARGET_CONTROL:
                acquisition = acquire_geo_series(accession)
                downloads.extend(acquisition["download_records"])
                _process_target_control(accession, acquisition, result)
        else:
            for accession in BENCHMARK:
                acquisition = acquire_geo_series(accession)
                downloads.extend(acquisition["download_records"])
                _process_benchmark(accession, acquisition, result)
    except StopAndReportError as exc:
        result.stop_conditions.append(exc.to_record())

    run_root = Path("runs") / "phase2" / track
    manifest = write_run_manifest(
        run_root,
        config=config,
        input_paths=[item["path"] for item in downloads],
        resource_records=downloads,
        package_names=PACKAGE_VERSION_NAMES + ["rpy2"],
        repo_root=".",
        allow_dirty=allow_dirty,
    )
    run_dirs = sorted(run_root.glob("*_*"))
    if run_dirs:
        summary_file = run_dirs[-1] / "phase2_result.json"
        summary_file.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n")
        result.files_created.append(str(summary_file))
        try:
            versions_file = run_dirs[-1] / "limma_versions.json"
            versions_file.write_text(json.dumps(limma_versions(), indent=2) + "\n")
            result.files_created.append(str(versions_file))
        except Exception:
            pass
    result.notes.append("Firewall: Phase 2 did not load disease-gene labels, drugs, or networks.")
    result.notes.append(f"plan_hash={manifest.get('plan_hash')}")
    _ = git
    _ = utc_now_path_component  # available if needed
    return result
