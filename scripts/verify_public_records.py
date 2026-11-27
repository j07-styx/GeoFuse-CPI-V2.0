#!/usr/bin/env python3
"""Check public V2.0 aggregate records, frozen assets and replay bundle identities.

This is a metadata/integrity check. It does not run inference or validate accuracy.
"""
from __future__ import annotations

import hashlib
import io
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


results = json.loads((ROOT / "results" / "PAPER_RESULTS.json").read_text(encoding="utf-8"))
native = json.loads((ROOT / "results" / "NATIVE_ACCURACY_AUDIT.json").read_text(encoding="utf-8"))
gate = json.loads((ROOT / "results" / "EXTERNAL_GATE.json").read_text(encoding="utf-8"))
identity = json.loads((ROOT / "manifests" / "frozen-artifact-identity.json").read_text(encoding="utf-8"))
run_manifest = json.loads((ROOT / "reproduction" / "replay-run-manifest.json").read_text(encoding="utf-8"))
openbind_report = json.loads((ROOT / "results" / "reproduction" / "openbind-reproduction.json").read_text(encoding="utf-8"))
biolip_report = json.loads((ROOT / "results" / "reproduction" / "biolip3-reproduction.json").read_text(encoding="utf-8"))

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

asset_path = ROOT / identity["inference_asset_archive"]
assert file_sha256(asset_path) == identity["inference_asset_archive_sha256"]
with zipfile.ZipFile(asset_path) as archive:
    names = set(archive.namelist())
    assert names == {"manifest.json", "GEOFUSE_CPI_V2.pt", "TRAIN_ANCHORS_NO_Y.npz"}
    assert "PRED_NOY.npz" not in names
    embedded = json.loads(archive.read("manifest.json"))
    for name, expected in embedded["files"].items():
        assert sha256(archive.read(name)) == expected

for name, record in identity["replay_bundles"].items():
    bundle_path = ROOT / record["path"]
    assert file_sha256(bundle_path) == record["sha256"], f"{name} replay bundle hash mismatch"
    with zipfile.ZipFile(bundle_path) as bundle:
        for member in bundle.namelist():
            if not member.endswith(".npz"):
                continue
            with zipfile.ZipFile(io.BytesIO(bundle.read(member))) as arrays:
                forbidden = {"y.npy", "label.npy", "labels.npy", "affinity.npy", "affinities.npy", "target.npy"}
                assert not (set(arrays.namelist()) & forbidden), f"Label-like array in {name}/{member}"

assert run_manifest["status"] == "PASS"
assert run_manifest["openbind"]["metric_checks"] == openbind_report["metric_checks"] == 60
assert run_manifest["openbind"]["bootstrap_interval_checks"] == openbind_report["bootstrap_interval_checks"] == 12
assert run_manifest["biolip3"]["metric_checks"] == biolip_report["metric_checks"] == 105
assert run_manifest["biolip3"]["bootstrap_interval_checks"] == biolip_report["bootstrap_interval_checks"] == 8
assert openbind_report["status"] == "PASS" and openbind_report["n"] == 482
assert biolip_report["status"] == "PASS"
assert identity["public_replay_verification"]["maximum_absolute_prediction_difference"] == 0.0

print("PASS: V2.0 aggregate metadata, frozen inference archive, replay bundle hashes and source-label replay reports are consistent.")
print("NOTE: this check does not rerun inference, source-label metrics or independent retraining.")
