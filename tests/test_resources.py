import pytest

from signflow.config import load_config
from signflow.resources import ResourceRecord, count_filtered_targets, verify_resources


def test_count_filtered_targets_uses_allowlist_and_threshold():
    rows = [
        {"datatypeScores": [{"id": "genetic_association", "score": 0.35}]},
        {"datatypeScores": [{"id": "genetic_association", "score": 0.19}]},
        {"datatypeScores": [{"id": "literature", "score": 0.99}]},
    ]
    assert count_filtered_targets(rows, {"genetic_association"}, 0.30) == 1


def test_resource_record_to_dict():
    record = ResourceRecord(
        name="demo",
        ok=True,
        status="ok",
        endpoint="https://example.org",
        checked_at_utc="2026-07-23T20:35:00Z",
    )
    assert record.to_dict()["name"] == "demo"


@pytest.mark.live
def test_verify_resources_live():
    config = load_config("configs/default.yaml")
    records = verify_resources(config)
    assert any(record["name"] == "open_targets" for record in records)
