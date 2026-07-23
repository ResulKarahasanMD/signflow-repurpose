from pathlib import Path

from signflow.config import load_config


def test_frozen_plan_constants_match_docs():
    config = load_config("configs/default.yaml")
    analysis_plan = Path("docs/phase0/ANALYSIS_PLAN.md").read_text()
    label_provenance = Path("docs/phase0/LABEL_PROVENANCE.md").read_text()
    adjudication = Path("docs/phase0/LABEL_ADJUDICATION.md").read_text()

    assert "`0.30`" in analysis_plan
    assert "**≥ 0.30**" in label_provenance
    assert "**N (minimum `plausible-ICBPS`)** | **15**" in adjudication
    assert "| `H7_FUTILITY_J` | `0.20` |" in analysis_plan
    assert "`MONDO_0005101`" in analysis_plan
    assert "`MONDO_0018301`" in analysis_plan
    assert "GSE75214`, `GSE87466" in analysis_plan
    assert "GSE11783`, `GSE57560" in analysis_plan
    assert "**396**" in label_provenance
    assert "**38**" in label_provenance

    plan = config["plan"]
    assert plan["label_score_threshold"] == 0.30
    assert plan["label_score_sensitivities"] == [0.20, 0.40]
    assert plan["conformal_alpha"] == 0.10
    assert plan["conformal_alpha_sensitivities"] == [0.05, 0.20]
    assert plan["null_permutations"] == 1000
    assert plan["signature_bootstraps"] == 1000
    assert plan["triangulation_k_primary"] == 20
    assert plan["triangulation_k_alternate"] == 50
    assert plan["random_seeds"] == [0, 1, 2, 3, 4]
    assert plan["module_sizes"] == [50, 100, 250]
    assert plan["label_adjudication_n"] == 15
    assert plan["h7_futility_j"] == 0.20
    assert plan["benchmark_disease_id"] == "MONDO_0005101"
    assert plan["target_disease_id"] == "MONDO_0018301"
    assert plan["benchmark_cohorts"] == ["GSE75214", "GSE87466"]
    assert plan["target_primary_cohorts"] == ["GSE11783", "GSE57560"]
    assert plan["target_inflammation_control"] == "GSE238208"
    assert plan["target_scrna"] == "GSE310557"
    assert plan["gse11783_primary_site"] == "non-ulcer"
    assert plan["uc_genetics_allowlist"] == ["genetic_association"]
    assert plan["uc_genetics_exclude"] == ["clinical", "literature", "rna_expression", "genetic_literature"]
    assert plan["benchmark_convergence_excludes"] == ["genetics"]
    assert plan["target_animal_model_status"] == "dropped"
    assert plan["frozen_counts"]["uc_genetics_primary"] == 396
    assert plan["frozen_counts"]["ic_animal_model_primary"] == 38
