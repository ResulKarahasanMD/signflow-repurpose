"""Configuration loading and runtime helpers for Phase 1."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import random
from pathlib import Path
from typing import Any

import numpy as np
import torch
import yaml


class ConfigError(ValueError):
    """Raised when configuration loading or overrides fail."""


def _merge_dicts(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_dicts(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _load_yaml_file(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text()) or {}
    if not isinstance(data, dict):
        raise ConfigError(f"Expected mapping at top level of {path}")
    return data


def _resolve_extends(parent: Path, extends: str | list[str]) -> list[Path]:
    items = [extends] if isinstance(extends, str) else list(extends)
    return [(parent / item).resolve() for item in items]


def _load_recursive(path: Path, seen: set[Path]) -> dict[str, Any]:
    resolved = path.resolve()
    if resolved in seen:
        raise ConfigError(f"Cyclic config extends detected at {resolved}")

    seen = set(seen)
    seen.add(resolved)
    data = _load_yaml_file(resolved)
    extends = data.pop("extends", None)

    if not extends:
        return data

    merged: dict[str, Any] = {}
    for base_path in _resolve_extends(resolved.parent, extends):
        merged = _merge_dicts(merged, _load_recursive(base_path, seen))
    return _merge_dicts(merged, data)


def dotted_get(data: dict[str, Any], key: str, default: Any | None = None) -> Any:
    current: Any = data
    for part in key.split("."):
        if not isinstance(current, dict) or part not in current:
            return default
        current = current[part]
    return current


def dotted_set(data: dict[str, Any], key: str, value: Any) -> None:
    current = data
    parts = key.split(".")
    for part in parts[:-1]:
        if part not in current or not isinstance(current[part], dict):
            current[part] = {}
        current = current[part]
    current[parts[-1]] = value


def parse_override(item: str) -> tuple[str, Any]:
    if "=" not in item:
        raise ConfigError(f"Override must be KEY=VALUE, got: {item}")
    key, raw_value = item.split("=", 1)
    key = key.strip()
    if not key:
        raise ConfigError(f"Override key cannot be empty: {item}")
    return key, yaml.safe_load(raw_value)


def apply_overrides(data: dict[str, Any], overrides: list[str] | None) -> dict[str, Any]:
    merged = copy.deepcopy(data)
    for item in overrides or []:
        key, value = parse_override(item)
        dotted_set(merged, key, value)
    return merged


def choose_device(cuda_available: bool, mps_available: bool) -> str:
    if cuda_available:
        return "cuda"
    if mps_available:
        return "mps"
    return "cpu"


def resolve_device(preferred: str = "auto") -> tuple[str, str | None]:
    preferred = preferred.lower()
    mps_available = bool(getattr(torch.backends, "mps", None)) and torch.backends.mps.is_available()

    if preferred == "auto":
        return choose_device(torch.cuda.is_available(), mps_available), None
    if preferred == "cuda":
        if torch.cuda.is_available():
            return "cuda", None
        return "cpu", "Requested device 'cuda' is unavailable; using cpu."
    if preferred == "mps":
        if mps_available:
            return "mps", None
        return "cpu", "Requested device 'mps' is unavailable; using cpu."
    if preferred == "cpu":
        return "cpu", None
    raise ConfigError(f"Unsupported device preference: {preferred}")


def get_seed_list(config: dict[str, Any]) -> list[int]:
    seeds = dotted_get(config, "plan.random_seeds")
    if not isinstance(seeds, list) or not all(isinstance(seed, int) for seed in seeds):
        raise ConfigError("plan.random_seeds must be a list of integers")
    return seeds


def seed_all(seed: int) -> dict[str, int]:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)

    seeded = {"python": seed, "python_hash_seed": seed, "numpy": seed, "torch": seed}
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        seeded["torch_cuda"] = seed
    if hasattr(torch.backends, "cudnn"):
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

    mps = getattr(torch, "mps", None)
    if mps is not None and hasattr(mps, "manual_seed"):
        mps.manual_seed(seed)
        seeded["torch_mps"] = seed

    return seeded


def canonical_config_json(config: dict[str, Any]) -> str:
    return json.dumps(config, sort_keys=True, separators=(",", ":"))


def semantic_config(config: dict[str, Any]) -> dict[str, Any]:
    semantic = copy.deepcopy(config)
    semantic.pop("runtime", None)
    return semantic


def plan_hash(config: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_config_json(semantic_config(config)).encode("utf-8")).hexdigest()


def resolved_config_hash(config: dict[str, Any]) -> str:
    return hashlib.sha256(canonical_config_json(config).encode("utf-8")).hexdigest()


def config_hash(config: dict[str, Any]) -> str:
    """Backward-compatible alias for the full resolved configuration hash."""
    return resolved_config_hash(config)


def load_config(config_path: str | Path, overrides: list[str] | None = None) -> dict[str, Any]:
    path = Path(config_path)
    if not path.exists():
        raise ConfigError(f"Config file does not exist: {path}")

    config = _load_recursive(path, seen=set())
    config = apply_overrides(config, overrides)
    runtime = config.setdefault("runtime", {})
    runtime.setdefault("python_version", "3.12")
    runtime.setdefault("device_preference", "auto")
    device, warning = resolve_device(runtime["device_preference"])
    runtime["device"] = device
    runtime["pytorch_enable_mps_fallback"] = os.environ.get("PYTORCH_ENABLE_MPS_FALLBACK", "1")
    runtime.setdefault(
        "pytorch_enable_mps_fallback_note",
        "Set PYTORCH_ENABLE_MPS_FALLBACK=1 when using Apple Silicon MPS.",
    )
    warnings = runtime.setdefault("warnings", [])
    if not isinstance(warnings, list):
        raise ConfigError("runtime.warnings must be a list")
    if warning:
        warnings.append(warning)
    return config
