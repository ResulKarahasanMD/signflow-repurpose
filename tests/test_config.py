from signflow.config import ConfigError, apply_overrides, choose_device, dotted_get, load_config, parse_override


def test_load_config_extends_and_defaults():
    config = load_config("configs/benchmark_uc.yaml")
    assert config["profile"]["track"] == "benchmark"
    assert config["plan"]["benchmark_cohorts"] == ["GSE75214", "GSE87466"]
    assert config["runtime"]["python_version"] == "3.12"
    assert config["runtime"]["device"] in {"cpu", "mps", "cuda"}


def test_dotted_get_and_overrides():
    config = load_config("configs/default.yaml")
    updated = apply_overrides(config, ["plan.label_score_threshold=0.4", "runtime.device_preference=cpu"])
    assert dotted_get(updated, "plan.label_score_threshold") == 0.4
    assert dotted_get(updated, "runtime.device_preference") == "cpu"


def test_parse_override_requires_key_value():
    try:
        parse_override("not-valid")
    except ConfigError:
        pass
    else:  # pragma: no cover
        raise AssertionError("Expected ConfigError")


def test_choose_device_priority():
    assert choose_device(cuda_available=True, mps_available=True) == "cuda"
    assert choose_device(cuda_available=False, mps_available=True) == "mps"
    assert choose_device(cuda_available=False, mps_available=False) == "cpu"
