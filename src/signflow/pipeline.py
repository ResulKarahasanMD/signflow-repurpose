"""Phase 1 pipeline entry points."""

from __future__ import annotations

from typing import Any

from .config import get_seed_list, seed_all


class StageNotImplementedError(RuntimeError):
    """Raised for intentionally deferred analysis stages."""


def run_smoke(config: dict[str, Any]) -> dict[str, Any]:
    seeds = get_seed_list(config)
    seeded = seed_all(seeds[0])
    return {
        "status": "ok",
        "device": config["runtime"]["device"],
        "seed_count": len(seeds),
        "first_seed": seeds[0],
        "seeded": seeded,
        "message": "Phase 1 smoke passed; analysis stages remain intentionally unimplemented.",
    }


def run_stage(stage_name: str, config: dict[str, Any]) -> dict[str, Any]:
    del config
    raise StageNotImplementedError(
        f"Stage '{stage_name}' is intentionally unimplemented in Phase 1. "
        "Phase 1 covers env/config/provenance/resource harness only."
    )
