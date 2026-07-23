"""Activity inference via decoupler CollecTRI and PROGENy."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import MANDATED_QC_HEADER, StopAndReportError, ensure_directory


def _score_with_ulm(expression: pd.DataFrame, network: pd.DataFrame) -> pd.DataFrame:
    import decoupler as dc

    scores = dc.mt.ulm(data=expression.T, net=network, tmin=5, verbose=False)
    # decoupler>=2 returns (estimate_df, pvalue_df)
    if isinstance(scores, tuple) and scores and isinstance(scores[0], pd.DataFrame):
        return scores[0]
    if isinstance(scores, pd.DataFrame):
        return scores
    if isinstance(scores, dict):
        for key in ("estimate", "score", "ulm_estimate"):
            if key in scores and isinstance(scores[key], pd.DataFrame):
                return scores[key]
        for value in scores.values():
            if isinstance(value, pd.DataFrame):
                return value
    raise StopAndReportError(
        f"decoupler returned an unexpected activity result object: {type(scores)}",
        reason="decoupler-resource-failure",
    )


def _dispersion_table(resource: str, activities: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "resource": resource,
            "source": activities.columns,
            "std": activities.std(axis=0).to_numpy(),
            "mad": (activities - activities.median(axis=0)).abs().median(axis=0).to_numpy(),
        }
    )


def infer_activities(expression: pd.DataFrame, output_dir: str | Path) -> dict[str, str]:
    """Infer TF/pathway activities and write dispersion QC tables."""

    try:
        import decoupler as dc

        collectri = dc.op.collectri(organism="human", remove_complexes=False, license="academic", verbose=False)
        progeny = dc.op.progeny(organism="human", top=100, license="academic", verbose=False)
        collectri_scores = _score_with_ulm(expression, collectri)
        progeny_scores = _score_with_ulm(expression, progeny)
    except Exception as exc:
        raise StopAndReportError(
            f"decoupler resource failure: {exc}",
            reason="decoupler-resource-failure",
        ) from exc

    output_dir = ensure_directory(output_dir)
    collectri_path = Path(output_dir) / "collectri_activities.tsv"
    progeny_path = Path(output_dir) / "progeny_activities.tsv"
    dispersion_path = Path(output_dir) / "activity_dispersion.tsv"
    header = f"# {MANDATED_QC_HEADER}\n"
    collectri_path.write_text(header)
    progeny_path.write_text(header)
    dispersion_path.write_text(header)
    collectri_scores.to_csv(collectri_path, sep="\t", mode="a")
    progeny_scores.to_csv(progeny_path, sep="\t", mode="a")
    pd.concat(
        [
            _dispersion_table("CollecTRI", collectri_scores),
            _dispersion_table("PROGENy", progeny_scores),
        ],
        ignore_index=True,
    ).to_csv(dispersion_path, sep="\t", mode="a", index=False)
    return {
        "collectri_activities": str(collectri_path),
        "progeny_activities": str(progeny_path),
        "activity_dispersion": str(dispersion_path),
    }
