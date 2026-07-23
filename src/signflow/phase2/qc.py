"""Descriptive QC outputs for Phase 2 expression studies."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA

from . import MANDATED_QC_HEADER, ensure_directory


def write_qc_outputs(
    expression: pd.DataFrame,
    phenotype_table: pd.DataFrame,
    study_id: str,
    output_dir: str | Path,
) -> list[str]:
    """Write lightweight descriptive QC outputs with the mandated header."""

    output_dir = ensure_directory(output_dir)
    files: list[str] = []

    header_path = Path(output_dir) / f"{study_id}_qc_header.txt"
    header_path.write_text(MANDATED_QC_HEADER + "\n")
    files.append(str(header_path))

    boxplot_path = Path(output_dir) / f"{study_id}_sample_boxplot.png"
    fig, ax = plt.subplots(figsize=(10, 4))
    expression.boxplot(ax=ax, rot=90)
    ax.set_title(MANDATED_QC_HEADER[:90] + "...")
    ax.set_ylabel("Expression")
    fig.tight_layout()
    fig.savefig(boxplot_path, dpi=150)
    plt.close(fig)
    files.append(str(boxplot_path))

    pca = PCA(n_components=2)
    scores = pca.fit_transform(expression.T.fillna(expression.T.mean()))
    pca_table = phenotype_table[["sample_id", "group"]].copy()
    pca_table["pc1"] = scores[:, 0]
    pca_table["pc2"] = scores[:, 1]
    pca_path = Path(output_dir) / f"{study_id}_pca.tsv"
    pca_table.to_csv(pca_path, sep="\t", index=False)
    files.append(str(pca_path))
    return files
