#!/usr/bin/env python3
"""Check consistency of the public aggregate-only V2.0 research records.

This script does not load model weights, perform inference, or validate accuracy.
"""
from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
results = json.loads((ROOT / "results" / "PAPER_RESULTS.json").read_text(encoding="utf-8"))
native = json.loads((ROOT / "results" / "NATIVE_ACCURACY_AUDIT.json").read_text(encoding="utf-8"))
gate = json.loads((ROOT / "results" / "EXTERNAL_GATE.json").read_text(encoding="utf-8"))
identity = json.loads((ROOT / "manifests" / "frozen-artifact-identity.json").read_text(encoding="utf-8"))

assert results["paper_model"] == "GeoFuse-CPI V2.0"
assert results["checkpoint_sha"] == identity["checkpoint_sha256"]
assert results["prediction_sha"] == identity["prediction_archive_sha256"]
assert native["checkpoint_sha"] == results["checkpoint_sha"]
assert native["prediction_sha"] == results["prediction_sha"]
assert gate["status"] == "BLOCKED_NO_ELIGIBLE_COHORT"
assert results["published_genSPARC_wins"] == 20

for key in ("davis__double_cold", "kiba__double_cold", "metz__double_cold"):
    row = results["metrics"][key]
    assert row["RMSE"] < row["matched_rows_reproduced_comparator"]["RMSE"]

print("PASS: public V2.0 metadata and aggregate-record identities are internally consistent.")
print("NOTE: this does not reproduce model predictions or independently validate benchmark accuracy.")
