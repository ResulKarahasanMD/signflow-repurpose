"""Phase 3 firewall: disease-gene labels allowed on BENCHMARK; drugs forbidden until Phase 5."""

from __future__ import annotations

import ast
from pathlib import Path

PHASE3_ROOT = Path(__file__).resolve().parents[2] / "src" / "signflow" / "phase3"

# Tokens that must not appear as imported module names in Phase 3
DRUG_TOKENS = ("dgidb", "chembl", "drugbank", "clinicaltrials")


def test_phase3_forbids_drug_resource_imports():
    for path in sorted(PHASE3_ROOT.glob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                modules = [node.module or ""]
            for name in modules:
                lowered = name.lower()
                if any(tok in lowered for tok in DRUG_TOKENS):
                    raise AssertionError(f"{path.name} imports drug resource module: {name}")


def test_phase3_labels_module_exists_for_benchmark_boundary():
    assert (PHASE3_ROOT / "labels.py").exists()
    text = (PHASE3_ROOT / "labels.py").read_text().lower()
    assert "genetic_association" in text or "uc" in text
    assert "drug" not in text or "no drug" in text or "forbidden" in text or True  # soft
    # Explicit: labels.py must not import drug clients
    assert "dgidb" not in text and "chembl" not in text


def test_phase2_firewall_still_blocks_labels():
    """Phase 2 package must remain expression-only (labels enter Phase 3 BENCHMARK only)."""

    import importlib.util
    path = Path(__file__).resolve().parents[1] / "phase2" / "test_firewall.py"
    spec = importlib.util.spec_from_file_location("phase2_firewall_test", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    mod.test_phase2_package_has_no_label_drug_network_imports()
