"""Phase 5 firewall: drugs allowed here; phase3 labels forbidden for scoring imports."""

from __future__ import annotations

import ast
from pathlib import Path

PHASE5 = Path(__file__).resolve().parents[2] / "src" / "signflow" / "phase5"
FORBIDDEN_LABEL_MODULES = ("signflow.phase3.labels",)


def test_phase5_does_not_import_phase3_labels():
    for path in PHASE5.glob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                mods = [node.module or ""]
            for m in mods:
                for bad in FORBIDDEN_LABEL_MODULES:
                    if m == bad or m.startswith(bad + "."):
                        raise AssertionError(f"{path.name} imports forbidden {m}")


def test_phase3_still_forbids_drug_modules():
    phase3 = Path(__file__).resolve().parents[2] / "src" / "signflow" / "phase3"
    drug_tokens = ("chembl", "dgidb", "drugbank")
    for path in phase3.glob("*.py"):
        text = path.read_text().lower()
        for tok in drug_tokens:
            if f"import {tok}" in text or f"from {tok}" in text:
                raise AssertionError(f"phase3/{path.name} mentions drug import {tok}")
