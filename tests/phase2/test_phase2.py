import gzip
from pathlib import Path

import pandas as pd
import pytest

from signflow.phase2 import StopAndReportError
from signflow.phase2.bootstrap_signature import stratified_bootstrap_samples
from signflow.phase2.meta_analysis import meta_analyze_signatures, random_effects_meta
from signflow.phase2.phenotype import classify_phenotypes, select_gse11783_primary_samples, select_uc_only
from signflow.phase2.process import collapse_probes_to_genes, detect_agilent_channel_count


def _write_gzip(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    return path


def test_ambiguous_phenotype_stops(tmp_path):
    matrix = _write_gzip(
        tmp_path / "GSE57560_series_matrix.txt.gz",
        """!Sample_title\t"s1"\n!Sample_geo_accession\t"GSM1"\n!Sample_characteristics_ch1\t"disease state: uncertain"\n!series_matrix_table_begin\n"ID_REF"\t"GSM1"\n"p1"\t1\n!series_matrix_table_end\n""",
    )
    with pytest.raises(StopAndReportError, match="could not be classified"):
        classify_phenotypes("GSE57560", matrix)


def test_mixed_ibd_uc_only_selection(tmp_path):
    matrix = _write_gzip(
        tmp_path / "GSE75214_series_matrix.txt.gz",
        (
            "!Sample_title\t\"uc\"\t\"ctrl\"\t\"cd\"\n"
            "!Sample_geo_accession\t\"GSM1\"\t\"GSM2\"\t\"GSM3\"\n"
            "!Sample_characteristics_ch1\t\"disease: ulcerative colitis\"\t\"disease: control\"\t\"disease: CD\"\n"
            "!Sample_characteristics_ch1\t\"disease activity: active\"\t\"disease activity: normal\"\t\"disease activity: active\"\n"
            "!Sample_characteristics_ch1\t\"tissue: colon\"\t\"tissue: colon\"\t\"tissue: ileum\"\n"
            "!series_matrix_table_begin\n"
            "\"ID_REF\"\t\"GSM1\"\t\"GSM2\"\t\"GSM3\"\n"
            "\"p1\"\t1\t2\t3\n"
            "!series_matrix_table_end\n"
        ),
    )
    phenotypes = classify_phenotypes("GSE75214", matrix)
    filtered, excluded = select_uc_only(phenotypes)
    assert set(filtered["group"]) == {"case", "control"}
    assert excluded["excluded"] == 1

def test_gse11783_primary_one_sample_per_patient(tmp_path):
    matrix = _write_gzip(
        tmp_path / "GSE11783_series_matrix.txt.gz",
        """!Sample_title\t"patient:12--location:ni--R:01"\t"patient:12--location:ulcus--R:01"\t"patient:13--location:ni--R:01"\t"patient:13--location:ulcus--R:01"\t"patient:2--location:healthy--R:01"\n!Sample_geo_accession\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\t"GSM5"\n!Sample_characteristics_ch1\t"disease state: IC, non-ulcer tissue"\t"disease state: IC, ulcer"\t"disease state: IC, non-ulcer tissue"\t"disease state: IC, ulcer"\t"disease state: no"\n!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\t"GSM3"\t"GSM4"\t"GSM5"\n"p1"\t1\t2\t3\t4\t5\n!series_matrix_table_end\n""",
    )
    phenotypes = classify_phenotypes("GSE11783", matrix)
    selected, mapping = select_gse11783_primary_samples(phenotypes)
    assert set(selected["sample_id"]) == {"GSM1", "GSM3", "GSM5"}
    assert mapping["selected_site"].tolist() == ["non-ulcer", "non-ulcer"]


def test_probe_collapse_rejects_non_protein_coding():
    expression = pd.DataFrame(
        {
            "probe_id": ["p1", "p2", "p3"],
            "GSM1": [5.0, 1.0, 7.0],
            "GSM2": [9.0, 1.1, 6.5],
        }
    )
    annotation = pd.DataFrame(
        {
            "ID": ["p1", "p2", "p3"],
            "Gene symbol": ["GENE1", "GENE1", "UOX"],
            "gene_type": ["protein_coding", "protein_coding", "lncRNA"],
        }
    )
    collapsed, qc = collapse_probes_to_genes(expression, annotation)
    assert list(collapsed.index) == ["GENE1"]
    assert collapsed.loc["GENE1", "GSM1"] == 5.0
    assert int(qc.loc[0, "rejected_non_protein_coding"]) == 1
    assert qc.loc[0, "collapse_rule"] == "max-mean probe per gene"


def test_meta_analysis_two_study_example():
    meta = random_effects_meta(pd.Series([1.0, 0.5]).to_numpy(), pd.Series([0.2, 0.3]).to_numpy())
    assert meta["n_studies"] == 2
    assert meta["q"] >= 0
    assert meta["i2"] >= 0

    combined = meta_analyze_signatures(
        [
            pd.DataFrame({"gene_symbol": ["GENE1"], "effect": [1.0], "se": [0.2], "study": ["A"]}),
            pd.DataFrame({"gene_symbol": ["GENE1"], "effect": [0.5], "se": [0.3], "study": ["B"]}),
        ]
    )
    assert combined.loc[0, "gene_symbol"] == "GENE1"
    assert combined.loc[0, "n_studies"] == 2


def test_stratified_bootstrap_preserves_group_sizes():
    phenotypes = pd.DataFrame(
        {
            "sample_id": ["case1", "case2", "ctrl1", "ctrl2"],
            "group": ["case", "case", "control", "control"],
        }
    )
    result = stratified_bootstrap_samples(phenotypes, n_bootstrap=25, seed=0)
    assert len(result["resamples"]) == 25
    for resample in result["resamples"]:
        assert len(resample) == 4
        assert sum(sample.startswith("case") for sample in resample) == 2
        assert sum(sample.startswith("ctrl") for sample in resample) == 2


def test_detect_agilent_channel_count(tmp_path):
    matrix = _write_gzip(
        tmp_path / "GSE57560_series_matrix.txt.gz",
        """!Sample_title\t"s1"\t"s2"\n!Sample_geo_accession\t"GSM1"\t"GSM2"\n!Sample_channel_count\t"1"\t"1"\n!Sample_characteristics_ch1\t"disease state: control"\t"disease state: Interstitial cystitis"\n!series_matrix_table_begin\n"ID_REF"\t"GSM1"\t"GSM2"\n"p1"\t1\t2\n!series_matrix_table_end\n""",
    )
    assert detect_agilent_channel_count(matrix) == 1


def test_i2_q_known_inputs():
    # Two identical effects -> Q=0, I2=0
    meta = random_effects_meta(__import__("numpy").asarray([1.0, 1.0]), __import__("numpy").asarray([0.2, 0.2]))
    assert meta["q"] == 0.0
    assert meta["i2"] == 0.0
    # Divergent effects -> positive Q/I2
    meta2 = random_effects_meta(__import__("numpy").asarray([1.0, -1.0]), __import__("numpy").asarray([0.1, 0.1]))
    assert meta2["q"] > 0
    assert meta2["i2"] > 0
