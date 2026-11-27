#!/usr/bin/env python3
"""Validate the aggregate ATOM3D historical record without opening source data.

This script checks the compact summary and provenance manifest only. It does not
read ATOM3D LMDB/parquet files, affinity labels, structures, feature arrays or
per-example predictions. It does not rerun or rescore the consumed official test.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
summary = json.loads((ROOT / "results" / "atom3d-historical-summary.json").read_text(encoding="utf-8"))
manifest = json.loads((ROOT / "manifests" / "atom3d-evidence-manifest.json").read_text(encoding="utf-8"))


def crosses_zero(interval: list[float]) -> bool:
    return len(interval) == 2 and interval[0] <= 0.0 <= interval[1]


assert summary["record_status"] == "HISTORICAL_ATOM3D_RESEARCH_NOT_V2_0_BENCHMARK"
assert summary["dataset"]["official_split_counts"] == {"train": 3507, "validation": 466, "test": 490}
assert summary["v2_0_boundary"]["atom3d_listed_in_v2_0_headline_ledger"] is False
assert summary["v2_0_boundary"]["frozen_v2_0_atom3d_evaluation_verified"] is False

validation = summary["official_validation_scored_once_by_earlier_model"]
assert validation["n"] == 466
assert validation["official_test_accessed_at_that_stage"] is False
assert crosses_zero(validation["native_minus_shuffled_group_bootstrap_ci95"])
assert crosses_zero(validation["native_minus_pocket_group_bootstrap_ci95"])

test = summary["official_test_scored_once_by_separate_earlier_model"]
assert test["n"] == 490
assert test["status"] == "ATOM3D_ORIGINAL_ID30_PHYSICS_TEST_SCORED_ONCE"
assert test["integrity_audit_no_test_reread_or_new_scoring"] is True
assert crosses_zero(test["native_minus_pocket_group_bootstrap_ci95"])
assert crosses_zero(test["native_minus_shuffled_group_bootstrap_ci95"])

assert manifest["official_test_integrity"]["no_test_reread_or_new_scoring"] is True
assert manifest["official_test_integrity"]["test_model_is_v2_0"] is False
assert manifest["redistribution_boundary"]["raw_dataset_included"] is False
assert manifest["redistribution_boundary"]["affinity_labels_included"] is False
assert manifest["redistribution_boundary"]["raw_structures_included"] is False
assert manifest["redistribution_boundary"]["derived_feature_arrays_included"] is False
assert manifest["redistribution_boundary"]["per_example_predictions_included"] is False
assert manifest["v2_0_headline_results_unchanged"] is True

branches = {item["branch"]: item for item in summary["historical_train_only_experiments"]}
assert {"V13", "V14", "V16", "V27", "V28_V29", "V30", "V31", "V33", "V34_V36"} <= set(branches)
assert branches["V31"]["native_minus_shuffled_protein_group_macro_rmse_ci95"][0] < 0
assert branches["V31"]["native_minus_shuffled_protein_group_macro_rmse_ci95"][1] > 0
assert branches["V33"]["strict_dual_cold_mse"]["native"] > branches["V33"]["strict_dual_cold_mse"]["shuffled"]

print("PASS: ATOM3D historical aggregate and provenance records are internally consistent.")
print("PASS: official test is marked consumed; no test data were opened or rescored by this verifier.")
print("NOTE: this is a record-integrity check, not a model evaluation or independent reproduction.")
