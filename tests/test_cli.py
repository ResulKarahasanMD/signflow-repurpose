from signflow.cli import main


def test_smoke_command_succeeds():
    assert main(["--allow-dirty", "smoke"]) == 0


def test_config_errors_exit_two():
    assert main(["--config", "missing.yaml", "smoke"]) == 2


def test_unimplemented_stage_exits_one():
    assert main(["--allow-dirty", "run", "--stage", "de"]) == 1
