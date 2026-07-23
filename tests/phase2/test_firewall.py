import ast
from pathlib import Path


PHASE2_ROOT = Path(__file__).resolve().parents[2] / "src" / "signflow" / "phase2"
ALLOWED_EXTERNAL = {
    "signflow.mapping",
    "signflow.provenance",
    "signflow.config",
    "signflow.resources",
}
DISALLOWED_TOKENS = ("label", "drug", "network")


def test_phase2_package_has_no_label_drug_network_imports():
    python_files = sorted(path for path in PHASE2_ROOT.glob("*.py") if path.name != "__init__.py")
    for path in python_files:
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                modules = [node.module or ""]
            for name in modules:
                lowered = name.lower()
                if any(token in lowered for token in DISALLOWED_TOKENS):
                    raise AssertionError(f"{path.name} violates the Phase 2 firewall: {name}")
                if name.startswith("signflow.") and not name.startswith("signflow.phase2"):
                    if name not in ALLOWED_EXTERNAL and not any(name.startswith(a + ".") for a in ALLOWED_EXTERNAL):
                        raise AssertionError(f"{path.name} imports disallowed external module: {name}")
