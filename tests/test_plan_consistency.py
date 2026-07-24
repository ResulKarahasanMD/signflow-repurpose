import re
from pathlib import Path

import yaml

from signflow.config import load_config

FROZEN_PLAN_PATH = Path("docs/phase0/FROZEN_PLAN.yaml")
ANALYSIS_PLAN_PATH = Path("docs/phase0/ANALYSIS_PLAN.md")
LABEL_PROVENANCE_PATH = Path("docs/phase0/LABEL_PROVENANCE.md")
LABEL_ADJUDICATION_PATH = Path("docs/phase0/LABEL_ADJUDICATION.md")


def _read_text(path: Path) -> str:
    return path.read_text()


def _read_frozen_plan() -> dict[str, object]:
    return yaml.safe_load(FROZEN_PLAN_PATH.read_text())["plan"]


def _require_match(pattern: str, text: str, key: str) -> re.Match[str]:
    match = re.search(pattern, text, re.MULTILINE)
    if not match:
        raise AssertionError(f"Missing doc regex match for {key}")
    return match


def _parse_backticked_list(cell: str) -> list[str]:
    values = re.findall(r"`([^`]+)`", cell)
    if values:
        return [value.strip().strip("{}") for value in values]
    values = [value.strip().strip("`").strip("{}") for value in cell.split(",") if value.strip()]
    if not values:
        raise AssertionError(f"Unable to parse list cell: {cell}")
    return values


def _assert_triplet(key: str, doc_value: object, frozen_plan: dict[str, object], config_plan: dict[str, object]) -> None:
    frozen_value = frozen_plan[key]
    config_value = config_plan[key]
    assert doc_value == frozen_value == config_value, (
        f"{key} mismatch: doc={doc_value!r} frozen={frozen_value!r} config={config_value!r}"
    )


def _assert_nested_triplet(
    parent_key: str,
    child_key: str,
    doc_value: object,
    frozen_plan: dict[str, object],
    config_plan: dict[str, object],
) -> None:
    frozen_value = frozen_plan[parent_key][child_key]
    config_value = config_plan[parent_key][child_key]
    assert doc_value == frozen_value == config_value, (
        f"{parent_key}.{child_key} mismatch: doc={doc_value!r} frozen={frozen_value!r} config={config_value!r}"
    )


def test_frozen_plan_matches_config_and_docs():
    frozen_plan = _read_frozen_plan()
    config_plan = load_config("configs/default.yaml")["plan"]
    analysis_plan = _read_text(ANALYSIS_PLAN_PATH)
    label_provenance = _read_text(LABEL_PROVENANCE_PATH)
    adjudication = _read_text(LABEL_ADJUDICATION_PATH)

    covered_keys: set[str] = set()
    for key, frozen_value in frozen_plan.items():
        assert key in config_plan, f"Config missing frozen key {key}"
        assert config_plan[key] == frozen_value, (
            f"{key} mismatch between frozen plan and config: "
            f"frozen={frozen_value!r} config={config_plan[key]!r}"
        )
        covered_keys.add(key)

    _assert_triplet(
        "label_score_threshold",
        float(_require_match(r"\| `LABEL_SCORE_THRESHOLD` \| `([0-9.]+)` \|", analysis_plan, "label_score_threshold").group(1)),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "label_score_sensitivities",
        [
            float(_require_match(r"\| Sensitivity low \| ≥ ([0-9.]+) \|", label_provenance, "label_score_sensitivities.low").group(1)),
            float(_require_match(r"\| Sensitivity high \| ≥ ([0-9.]+) \|", label_provenance, "label_score_sensitivities.high").group(1)),
        ],
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "random_seeds",
        [
            int(value)
            for value in _require_match(r"\| `RANDOM_SEEDS` \| `\{([0-9,]+)\}` \|", analysis_plan, "random_seeds")
            .group(1)
            .split(",")
        ],
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "label_adjudication_n",
        int(
            _require_match(
                r"\| \*\*N \(minimum `plausible-ICBPS`\)\*\* \| \*\*([0-9]+)\*\* \|",
                adjudication,
                "label_adjudication_n",
            ).group(1)
        ),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "h7_futility_j",
        float(_require_match(r"\| `H7_FUTILITY_J` \| `([0-9.]+)` \|", analysis_plan, "h7_futility_j").group(1)),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "benchmark_disease_id",
        _require_match(r"\| `BENCHMARK_DISEASE_ID` \| `([^`]+)` \|", analysis_plan, "benchmark_disease_id").group(1),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "target_disease_id",
        _require_match(r"\| `TARGET_DISEASE_ID` \| `([^`]+)` \|", analysis_plan, "target_disease_id").group(1),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "benchmark_cohorts",
        _parse_backticked_list(
            _require_match(r"\| `BENCHMARK_COHORTS` \| ([^|]+) \|", analysis_plan, "benchmark_cohorts").group(1)
        ),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "target_primary_cohorts",
        _parse_backticked_list(
            _require_match(r"\| `TARGET_PRIMARY_COHORTS` \| ([^|]+) \|", analysis_plan, "target_primary_cohorts").group(1)
        ),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "target_inflammation_control",
        _require_match(r"\| `TARGET_INFLAMMATION_CONTROL` \| `([^`]+)` \|", analysis_plan, "target_inflammation_control").group(1),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "target_scrna",
        _require_match(r"\| `TARGET_SCRNA` \| `([^`]+)` \|", analysis_plan, "target_scrna").group(1),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "gse11783_primary_site",
        _require_match(r"\| `GSE11783_PRIMARY_SITE` \| `([^`]+)` \|", analysis_plan, "gse11783_primary_site").group(1),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "uc_genetics_allowlist",
        _parse_backticked_list(
            _require_match(
                r"\| BENCHMARK \| `MONDO_0005101` \| ([^|]+) \| [^|]+ \|",
                label_provenance,
                "uc_genetics_allowlist",
            ).group(1)
        ),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "uc_genetics_exclude",
        _parse_backticked_list(
            _require_match(
                r"\| BENCHMARK \| `MONDO_0005101` \| [^|]+ \| ([^|]+) \|",
                label_provenance,
                "uc_genetics_exclude",
            ).group(1)
        ),
        frozen_plan,
        config_plan,
    )
    _assert_triplet(
        "target_animal_model_status",
        "dropped"
        if _require_match(
            r"DROP `animal_model` from TARGET labels entirely",
            adjudication,
            "target_animal_model_status",
        )
        else "missing",
        frozen_plan,
        config_plan,
    )
    _assert_nested_triplet(
        "frozen_counts",
        "uc_genetics_primary",
        int(
            _require_match(
                r"\| UC \| `genetic_association` only \| [0-9]+ \| \*\*([0-9]+)\*\* \| [0-9]+ \|",
                label_provenance,
                "frozen_counts.uc_genetics_primary",
            ).group(1)
        ),
        frozen_plan,
        config_plan,
    )
    _assert_nested_triplet(
        "frozen_counts",
        "ic_animal_model_primary",
        int(
            _require_match(
                r"\| IC/BPS \| `animal_model` only \| [0-9]+ \| \*\*([0-9]+)\*\* \| [0-9]+ \|",
                label_provenance,
                "frozen_counts.ic_animal_model_primary",
            ).group(1)
        ),
        frozen_plan,
        config_plan,
    )

    assert covered_keys == set(frozen_plan), (
        f"Frozen plan keys missing assertion coverage: {sorted(set(frozen_plan) - covered_keys)}"
    )

    # Gold-standard drug file hash must match frozen digest (Phase 5 §B)
    import hashlib
    gold_path = Path(frozen_plan["benchmark_gold_drugs_path"])
    assert gold_path.exists(), f"Missing gold drugs file {gold_path}"
    digest = hashlib.sha256(gold_path.read_bytes()).hexdigest()
    assert digest == frozen_plan["benchmark_gold_drugs_sha256"], (
        f"benchmark_gold_drugs_sha256 mismatch: file={digest} frozen={frozen_plan['benchmark_gold_drugs_sha256']}"
    )
    assert frozen_plan["benchmark_gold_out_of_substrate_policy"] == (
        "exclude_from_denominator_report_as_coverage_miss"
    )

    # Phase 5.3 / 5.4 freezes
    assert frozen_plan["phase5_score_variant_primary"] == "signed"
    assert frozen_plan["phase5_score_variant_coprimary"] == "unsigned"
    assert frozen_plan["phase5_score_variant_unsigned_matched"] == "unsigned_matched"
    assert frozen_plan["phase5_score_variant_unsigned_full"] == "unsigned_full"
    assert frozen_plan["phase5_h7_top_k"] == 20
    assert frozen_plan["phase5_h7_n_bootstrap"] == 1000
    assert frozen_plan["phase5_control_graph_variants"] == ["rewired", "edge_shuffled", "sign_shuffled"]
    assert Path(frozen_plan["phase5_soft_positives_path"]).exists()
    assert frozen_plan["phase5_candidate_pool_rule"] == "disease_agnostic_approved_drugs"
    assert frozen_plan["phase5_candidate_pool_disease_conditioning"] == "forbidden"
    assert frozen_plan["phase5_candidate_pool_min_in_substrate"] == 500
    assert frozen_plan["phase5_candidate_pool_gold_prevalence_ceiling"] == 0.10
    assert "ion_channel" in frozen_plan["phase5_target_class_taxonomy"]
    assert frozen_plan["phase5_score_k_values"] == [20, 50, 100, 250]
    assert frozen_plan["phase5_candidate_pool_sha256"] == (
        "f748e4b93a24bd491d6040b0797ed081adc9b010484d19532ff3b4e4261d294c"
    )
    assert frozen_plan["phase5_candidate_pool_n_universe"] == 4225
    assert frozen_plan["phase5_candidate_pool_n_gold"] == 29
    pool_path = Path(frozen_plan["phase5_candidate_pool_path"])
    if pool_path.exists():
        import hashlib

        digest = hashlib.sha256(pool_path.read_bytes()).hexdigest()
        assert digest == frozen_plan["phase5_candidate_pool_sha256"]
