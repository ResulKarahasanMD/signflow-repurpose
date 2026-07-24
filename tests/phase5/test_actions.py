"""Phase 5 action vocabulary tests."""

from __future__ import annotations

from signflow.phase5.actions import UNKNOWN, action_sign, is_signed_action, normalize_action_token


def test_normalize_inhibitor_variants():
    assert normalize_action_token("Inhibition") == "inhibitor"
    assert normalize_action_token("negative modulator") == "inhibitor"


def test_unknown_action_not_signed():
    assert normalize_action_token(None) == UNKNOWN
    assert normalize_action_token("binder") == UNKNOWN
    assert action_sign("binder") is None
    assert not is_signed_action("binder")


def test_signed_actions_map_to_signs():
    assert action_sign("agonist") == 1
    assert action_sign("inhibitor") == -1
