import copy

import pytest

from signflow.config import (
    ConfigError,
    apply_overrides,
    choose_device,
    dotted_get,
    load_config,
    parse_override,
    plan_hash,
    resolve_device,
    resolved_config_hash,
    semantic_config,
)


def test_load_config_extends_and_defaults():
    config = load_config("configs/benchmark_uc.yaml")
    assert config["profile"]["track"] == "benchmark"
    assert config["plan"]["benchmark_cohorts"] == ["GSE75214", "GSE87466"]
    assert config["plan"]["triangulation_excludes"] == ["genetics"]
    assert config["runtime"]["python_version"] == "3.12"
    assert config["runtime"]["device"] in {"cpu", "mps", "cuda"}
    assert config["runtime"]["warnings"] == []


def test_dotted_get_and_overrides():
    config = load_config("configs/default.yaml")
    updated = apply_overrides(config, ["plan.label_score_threshold=0.4", "runtime.device_preference=cpu"])
    assert dotted_get(updated, "plan.label_score_threshold") == 0.4
    assert dotted_get(updated, "runtime.device_preference") == "cpu"


def test_parse_override_requires_key_value():
    with pytest.raises(ConfigError):
        parse_override("not-valid")


def test_choose_device_priority():
    assert choose_device(cuda_available=True, mps_available=True) == "cuda"
    assert choose_device(cuda_available=False, mps_available=True) == "mps"
    assert choose_device(cuda_available=False, mps_available=False) == "cpu"


def test_resolve_device_warns_when_requested_accelerator_unavailable(monkeypatch):
    monkeypatch.setattr("signflow.config.torch.cuda.is_available", lambda: False)
    monkeypatch.setattr("signflow.config.torch.backends.mps.is_available", lambda: False)
    device, warning = resolve_device("cuda")
    assert device == "cpu"
    assert warning == "Requested device 'cuda' is unavailable; using cpu."


def test_plan_hash_excludes_runtime_changes():
    config = load_config("configs/default.yaml")
    mutated = copy.deepcopy(config)
    mutated["runtime"]["device"] = "cuda"
    mutated["runtime"]["warnings"].append("synthetic warning")

    assert semantic_config(config) == semantic_config(mutated)
    assert plan_hash(config) == plan_hash(mutated)
    assert resolved_config_hash(config) != resolved_config_hash(mutated)
