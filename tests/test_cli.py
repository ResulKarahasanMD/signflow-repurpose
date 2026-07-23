import subprocess
from pathlib import Path

from signflow.cli import main


def _init_git_repo(repo: Path) -> None:
    subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "config", "user.email", "tests@example.com"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "config", "user.name", "Test User"], cwd=repo, check=True, capture_output=True, text=True)
    (repo / "README.md").write_text("fixture\n")
    subprocess.run(["git", "add", "README.md"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "commit", "-m", "init"], cwd=repo, check=True, capture_output=True, text=True)


def test_smoke_command_succeeds(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    repo.mkdir()
    _init_git_repo(repo)
    monkeypatch.chdir(repo)

    config_path = (Path(__file__).resolve().parents[1] / "configs" / "default.yaml").as_posix()
    assert main(["--config", config_path, "smoke"]) == 0
    manifests = list((repo / "runs" / "smoke").glob("*/run_manifest.json"))
    assert len(manifests) == 1


def test_config_errors_exit_two():
    assert main(["--config", "missing.yaml", "smoke"]) == 2


def test_unimplemented_stage_exits_one():
    config_path = (Path(__file__).resolve().parents[1] / "configs" / "default.yaml").as_posix()
    assert main(["--config", config_path, "run", "--stage", "de"]) == 1
