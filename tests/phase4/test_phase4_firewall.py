"""Phase 4 firewall: no drug imports; labels BENCHMARK-only path."""

from __future__ import annotations

import ast
from pathlib import Path

PHASE4 = Path(__file__).resolve().parents[2] / "src" / "signflow" / "phase4"
DRUG = ("dgidb", "chembl", "drugbank", "clinicaltrials")


def test_phase4_forbids_drug_imports():
    for path in PHASE4.glob("*.py"):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                mods = [node.module or ""]
            for m in mods:
                if any(d in m.lower() for d in DRUG):
                    raise AssertionError(f"{path.name} imports drug module {m}")
