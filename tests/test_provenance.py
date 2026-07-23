import subprocess

import pytest

from signflow.provenance import DirtyRepositoryError, GitStateError, git_state, write_run_manifest


def _init_git_repo(repo):
    subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "config", "user.email", "tests@example.com"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "config", "user.name", "Test User"], cwd=repo, check=True, capture_output=True, text=True)
    (repo / "README.md").write_text("fixture\n")
    subprocess.run(["git", "add", "README.md"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "commit", "-m", "init"], cwd=repo, check=True, capture_output=True, text=True)


def _config():
    return {
        "plan": {"random_seeds": [0, 1, 2, 3, 4]},
        "runtime": {"device": "cpu", "warnings": []},
    }


def test_git_state_raises_outside_git_repo(tmp_path):
    with pytest.raises(GitStateError):
        git_state(tmp_path)


def test_write_run_manifest_refuses_dirty_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _init_git_repo(repo)
    (repo / "input.txt").write_text("hello")

    with pytest.raises(DirtyRepositoryError):
        write_run_manifest(
            repo / "runs" / "smoke",
            config=_config(),
            input_paths=[repo / "input.txt"],
            resource_records=[],
            package_names=["numpy"],
            repo_root=repo,
            allow_dirty=False,
        )


def test_write_run_manifest_allows_dirty_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _init_git_repo(repo)
    input_path = repo / "input.txt"
    input_path.write_text("hello")

    manifest = write_run_manifest(
        repo / "runs" / "smoke",
        config=_config(),
        input_paths=[input_path],
        resource_records=[{"name": "demo", "ok": True}],
        package_names=["numpy"],
        repo_root=repo,
        allow_dirty=True,
    )
    assert manifest["git"]["dirty"] is True
    assert manifest["runtime_environment"]["device"] == "cpu"
    assert manifest["input_sha256"][str(input_path)]


def test_write_run_manifest_is_append_only(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _init_git_repo(repo)
    input_path = repo / "input.txt"
    input_path.write_text("hello")
    subprocess.run(["git", "add", "input.txt"], cwd=repo, check=True, capture_output=True, text=True)
    subprocess.run(["git", "commit", "-m", "add input"], cwd=repo, check=True, capture_output=True, text=True)

    output_root = repo / "runs" / "smoke"
    write_run_manifest(
        output_root,
        config=_config(),
        input_paths=[input_path],
        resource_records=[],
        package_names=["numpy"],
        repo_root=repo,
        allow_dirty=True,
    )
    write_run_manifest(
        output_root,
        config=_config(),
        input_paths=[input_path],
        resource_records=[],
        package_names=["numpy"],
        repo_root=repo,
        allow_dirty=True,
    )

    manifest_dirs = sorted(path for path in output_root.iterdir() if path.is_dir())
    assert len(manifest_dirs) == 2
    assert all((path / "run_manifest.json").exists() for path in manifest_dirs)
