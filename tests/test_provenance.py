import subprocess

from signflow.provenance import DirtyRepositoryError, write_run_manifest


def _init_git_repo(repo):
    subprocess.run(["git", "init"], cwd=repo, check=True, capture_output=True, text=True)


def test_write_run_manifest_refuses_dirty_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _init_git_repo(repo)
    (repo / "input.txt").write_text("hello")

    config = {"plan": {"random_seeds": [0, 1, 2, 3, 4]}}
    try:
        write_run_manifest(
            repo / "run_manifest.json",
            config=config,
            input_paths=[repo / "input.txt"],
            resource_records=[],
            package_names=["numpy"],
            repo_root=repo,
            allow_dirty=False,
        )
    except DirtyRepositoryError:
        pass
    else:  # pragma: no cover
        raise AssertionError("Expected dirty repository refusal")


def test_write_run_manifest_allows_dirty_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _init_git_repo(repo)
    input_path = repo / "input.txt"
    input_path.write_text("hello")

    config = {"plan": {"random_seeds": [0, 1, 2, 3, 4]}}
    manifest = write_run_manifest(
        repo / "run_manifest.json",
        config=config,
        input_paths=[input_path],
        resource_records=[{"name": "demo", "ok": True}],
        package_names=["numpy"],
        repo_root=repo,
        allow_dirty=True,
    )
    assert manifest["git"]["dirty"] is True
    assert manifest["input_sha256"][str(input_path)]
