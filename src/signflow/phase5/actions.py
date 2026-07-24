"""Normalize drug action strings via frozen phase5_action_sign_map."""

from __future__ import annotations

import yaml

from . import FROZEN_PLAN

UNKNOWN = "unknown_action"


def load_action_sign_map() -> dict[str, int | None]:
    plan = yaml.safe_load(FROZEN_PLAN.read_text())["plan"]
    raw = plan["phase5_action_sign_map"]
    return {str(k): (None if v is None else int(v)) for k, v in raw.items()}


def normalize_action_token(raw: str | None) -> str:
    """Map free-text / enum action to canonical action key."""

    if raw is None:
        return UNKNOWN
    text = str(raw).strip().lower()
    if not text or text in {"nan", "none", "unknown", "other/unknown"}:
        return UNKNOWN
    mapping = {
        "inhibitor": "inhibitor",
        "inhibition": "inhibitor",
        "inhibitory": "inhibitor",
        "antagonist": "antagonist",
        "antagonism": "antagonist",
        "agonist": "agonist",
        "activation": "activator",
        "activator": "activator",
        "positive modulator": "activator",
        "negative modulator": "inhibitor",
        "blocker": "inhibitor",
        "binder": UNKNOWN,
        "modulator": UNKNOWN,
    }
    if text in mapping:
        return mapping[text]
    for key in ("agonist", "activator", "antagonist", "inhibitor"):
        if key in text:
            return key
    return UNKNOWN


def action_sign(raw: str | None, sign_map: dict[str, int | None] | None = None) -> int | None:
    sign_map = sign_map or load_action_sign_map()
    canonical = normalize_action_token(raw)
    return sign_map.get(canonical, sign_map.get(UNKNOWN))


def is_signed_action(raw: str | None, sign_map: dict[str, int | None] | None = None) -> bool:
    return action_sign(raw, sign_map=sign_map) is not None
