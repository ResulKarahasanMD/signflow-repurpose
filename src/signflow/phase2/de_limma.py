"""Differential expression via rpy2 and Bioconductor limma."""

from __future__ import annotations

import importlib.util
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from . import StopAndReportError

LIMMA_INSTALL_HINT = (
    "R + Bioconductor limma are required for Phase 2 DE. Install rpy2, ensure R/Rscript are on PATH, "
    "and run: Rscript -e 'if (!requireNamespace(\"BiocManager\", quietly=TRUE)) "
    "install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org\"); "
    "BiocManager::install(\"limma\", ask=FALSE, update=FALSE)'"
)


def _require_r_limma() -> None:
    if not importlib.util.find_spec("rpy2"):
        raise StopAndReportError(f"rpy2 is not installed. {LIMMA_INSTALL_HINT}", reason="r-limma-missing")

    command = ["Rscript", "-e", 'quit(status=ifelse(requireNamespace("limma", quietly=TRUE), 0, 1))']
    probe = subprocess.run(command, capture_output=True, text=True, check=False)
    if probe.returncode != 0:
        raise StopAndReportError(f"Bioconductor limma is not installed. {LIMMA_INSTALL_HINT}", reason="r-limma-missing")


def standardized_effects(expression: pd.DataFrame, phenotype_table: pd.DataFrame) -> pd.DataFrame:
    """Unmoderated Hedges' g — sensitivity only (not primary meta input)."""

    case_ids = phenotype_table.loc[phenotype_table["group"] == "case", "sample_id"].tolist()
    control_ids = phenotype_table.loc[phenotype_table["group"] == "control", "sample_id"].tolist()
    case_values = expression[case_ids].to_numpy(dtype=float)
    control_values = expression[control_ids].to_numpy(dtype=float)
    n_case = len(case_ids)
    n_control = len(control_ids)
    case_mean = np.nanmean(case_values, axis=1)
    control_mean = np.nanmean(control_values, axis=1)
    case_var = np.nanvar(case_values, axis=1, ddof=1)
    control_var = np.nanvar(control_values, axis=1, ddof=1)
    pooled_sd = np.sqrt(((n_case - 1) * case_var + (n_control - 1) * control_var) / max(n_case + n_control - 2, 1))
    pooled_sd = np.where(pooled_sd == 0, np.nan, pooled_sd)
    cohen_d = (case_mean - control_mean) / pooled_sd
    correction = 1.0 - (3.0 / max(4 * (n_case + n_control) - 9, 1))
    hedges_g = cohen_d * correction
    se = np.sqrt(((n_case + n_control) / (n_case * n_control)) + (hedges_g**2 / max(2 * (n_case + n_control - 2), 1)))
    return pd.DataFrame(
        {
            "gene_symbol": expression.index,
            "hedges_g": hedges_g,
            "hedges_g_se": se,
            "n_case": n_case,
            "n_control": n_control,
        }
    )


def limma_significance_counts(de_table: pd.DataFrame) -> dict[str, int]:
    fdr = de_table["adj.P.Val"] if "adj.P.Val" in de_table.columns else de_table.get("bh_fdr")
    if fdr is None:
        return {"bh_fdr_lt_0.05": 0, "bh_fdr_lt_0.10": 0, "n_genes": int(len(de_table))}
    return {
        "bh_fdr_lt_0.05": int((fdr < 0.05).sum()),
        "bh_fdr_lt_0.10": int((fdr < 0.10).sum()),
        "n_genes": int(len(de_table)),
        "nan_moderated_effect": int(de_table["effect"].isna().sum()) if "effect" in de_table.columns else 0,
        "nan_moderated_se": int(de_table["se"].isna().sum()) if "se" in de_table.columns else 0,
    }


def run_limma_differential_expression(
    expression: pd.DataFrame,
    phenotype_table: pd.DataFrame,
    *,
    accession: str,
) -> pd.DataFrame:
    """Run limma moderated DE; primary effect/se are moderated logFC and moderated SE."""

    _require_r_limma()

    required_columns = {"sample_id", "group"}
    if not required_columns.issubset(phenotype_table.columns):
        raise StopAndReportError(
            f"Phenotype table for {accession} is missing required DE columns",
            accession=accession,
            reason="ambiguous-phenotype",
        )

    sample_ids = phenotype_table["sample_id"].tolist()
    if set(phenotype_table["group"]) != {"case", "control"}:
        raise StopAndReportError(
            f"{accession} DE requires exactly case/control samples after filtering",
            accession=accession,
            reason="ambiguous-phenotype",
        )

    expr = expression.loc[:, sample_ids]
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
    except Exception as exc:  # pragma: no cover
        raise StopAndReportError(
            f"rpy2 import failed after dependency check: {exc}. {LIMMA_INSTALL_HINT}",
            reason="r-limma-missing",
        ) from exc

    expr_df = expr.copy()
    expr_df.insert(0, "gene_symbol", expr_df.index)
    sample_df = phenotype_table[["sample_id", "group"]].copy()
    with localconverter(ro.default_converter + pandas2ri.converter):
        ro.globalenv["expr_df"] = ro.conversion.get_conversion().py2rpy(expr_df)
        ro.globalenv["sample_df"] = ro.conversion.get_conversion().py2rpy(sample_df)
        result = ro.r(
            """
            suppressPackageStartupMessages(library(limma))
            rownames(expr_df) <- expr_df$gene_symbol
            expr_df$gene_symbol <- NULL
            expr_mat <- as.matrix(expr_df)
            sample_df$group <- factor(sample_df$group, levels=c("control", "case"))
            rownames(sample_df) <- sample_df$sample_id
            design <- model.matrix(~ 0 + group, data=sample_df)
            colnames(design) <- c("control", "case")
            contrast <- makeContrasts(case - control, levels=design)
            fit <- lmFit(expr_mat, design)
            fit <- contrasts.fit(fit, contrast)
            fit <- eBayes(fit, robust=TRUE)
            tt <- topTable(fit, number=Inf, sort.by="none")
            coef <- 1
            tt$moderated_logFC <- as.numeric(fit$coefficients[, coef])
            tt$moderated_se <- as.numeric(fit$stdev.unscaled[, coef] * sqrt(fit$s2.post))
            tt$moderated_t <- as.numeric(fit$t[, coef])
            tt$df_total <- as.numeric(fit$df.total)
            tt$gene_symbol <- rownames(tt)
            tt
            """
        )
        de_table = ro.conversion.get_conversion().rpy2py(result)

    # Primary meta inputs: moderated logFC + moderated SE
    de_table = de_table.rename(columns={})
    de_table["effect"] = de_table["moderated_logFC"]
    de_table["se"] = de_table["moderated_se"]
    # Sensitivity: unmoderated Hedges' g
    g_table = standardized_effects(expr, phenotype_table)
    merged = de_table.merge(g_table, on="gene_symbol", how="left")
    # Drop genes with non-finite moderated SE before meta (counted in QC)
    return merged.sort_values("P.Value").reset_index(drop=True)


def limma_versions() -> dict[str, str]:
    _require_r_limma()
    import rpy2.robjects as ro

    return {
        "r_version": str(ro.r("R.version.string")[0]),
        "limma_version": str(ro.r('as.character(packageVersion("limma"))')[0]),
    }


def write_de_plots(de_table: pd.DataFrame, output_dir: str | Path, *, accession: str) -> list[str]:
    import matplotlib.pyplot as plt

    from . import ensure_directory

    output_dir = ensure_directory(output_dir)
    files: list[str] = []
    logfc = de_table["logFC"] if "logFC" in de_table.columns else de_table["effect"]
    pvals = de_table["P.Value"] if "P.Value" in de_table.columns else pd.Series(np.nan, index=de_table.index)
    neglog = -np.log10(pvals.clip(lower=1e-300))

    volcano_path = Path(output_dir) / f"{accession}_volcano.png"
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(logfc, neglog, s=6, alpha=0.4, c="0.3")
    ax.set_xlabel("moderated log2FC")
    ax.set_ylabel("-log10(P)")
    ax.set_title(f"{accession} volcano")
    fig.tight_layout()
    fig.savefig(volcano_path, dpi=150)
    plt.close(fig)
    files.append(str(volcano_path))

    if "AveExpr" in de_table.columns:
        ma_path = Path(output_dir) / f"{accession}_ma.png"
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(de_table["AveExpr"], logfc, s=6, alpha=0.4, c="0.3")
        ax.axhline(0, color="0.5", lw=1)
        ax.set_xlabel("AveExpr")
        ax.set_ylabel("moderated log2FC")
        ax.set_title(f"{accession} MA")
        fig.tight_layout()
        fig.savefig(ma_path, dpi=150)
        plt.close(fig)
        files.append(str(ma_path))
    return files
