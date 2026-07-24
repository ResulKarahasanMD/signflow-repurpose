"""drug_data_touched must be true whenever Phase 5 reads drug resources."""

from __future__ import annotations

from pathlib import Path

import pytest

from signflow.phase5 import StopAndReportError
from signflow.phase5.run import _mark_drug_data_touched, assert_drug_data_touched

PHASE5 = Path(__file__).resolve().parents[2] / "src" / "signflow" / "phase5"


def test_mark_drug_data_touched_sets_true_when_resources_read():
    result: dict[str, object] = {}
    _mark_drug_data_touched(result, resources_read=True)
    assert result["drug_data_touched"] is True


def test_assert_drug_data_touched_raises_when_false_after_read():
    with pytest.raises(StopAndReportError) as exc:
        assert_drug_data_touched({"drug_data_touched": False}, resources_read=True)
    assert exc.value.reason == "drug-data-touched-false"


def test_assert_drug_data_touched_ok_when_true():
    assert_drug_data_touched({"drug_data_touched": True}, resources_read=True)


def test_run_phase5_source_forces_drug_data_touched_true():
    text = (PHASE5 / "run.py").read_text()
    assert "assert_drug_data_touched" in text
    assert "resources_read=True" in text


def test_acquire_marks_drug_data_touched():
    text = (PHASE5 / "acquire.py").read_text()
    assert '"drug_data_touched": True' in text
