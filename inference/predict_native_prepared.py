#!/usr/bin/env python3
"""Native-expert inference from already prepared, label-free 423D/704D features."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from geofuse_native_expert import predict_native

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = ROOT / "models" / "v2.0" / "GEOFUSE_CPI_V2.pt"
CHECKPOINT_SHA256 = "b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, type=Path,
                        help="Label-free NPZ with x2d (N,423), real/shuffled/pocket (N,704), optional ids.")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT)
    args = parser.parse_args()

    if args.out.exists():
        parser.error(f"Refusing to overwrite existing output: {args.out}")
    if not args.features.is_file() or not args.checkpoint.is_file():
        parser.error("Input feature archive or frozen checkpoint is missing.")
    if sha256_file(args.checkpoint) != CHECKPOINT_SHA256:
        parser.error("Checkpoint SHA-256 does not match the frozen V2.0 release.")

    with np.load(args.features, allow_pickle=False) as archive:
        if any(key.lower() in {"y", "label", "labels", "affinity", "affinities", "target"} for key in archive.files):
            raise ValueError("Inference archive must be label-free; affinity/target fields are forbidden.")
        required = {"x2d", "real", "shuffled", "pocket"}
        if not required.issubset(archive.files):
            raise ValueError(f"Missing feature arrays: {sorted(required - set(archive.files))}")
        x2d = archive["x2d"].astype(np.float32)
        geometry = {key: archive[key].astype(np.float32) for key in ("real", "shuffled", "pocket")}
        ids = archive["ids"].astype(str) if "ids" in archive.files else np.arange(len(x2d)).astype(str)

    n = len(x2d)
    if x2d.shape != (n, 423) or any(value.shape != (n, 704) for value in geometry.values()):
        raise ValueError("Expected x2d=(N,423) and each geometry array=(N,704).")
    if len(ids) != n or len(set(ids.tolist())) != n:
        raise ValueError("IDs must be unique and aligned with the feature rows.")
    if not np.isfinite(x2d).all() or not all(np.isfinite(value).all() for value in geometry.values()):
        raise ValueError("Nonfinite prepared feature values.")

    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=True)
    outputs = predict_native(checkpoint, x2d, geometry)
    if any(not np.isfinite(np.asarray(outputs[key])).all() for key in ("2d", "real", "shuffled", "pocket", "raw")):
        raise RuntimeError("Native expert returned nonfinite predictions.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["sample_id", "calibrated_2d", "native_real", "shuffled_geometry_control", "pocket_only_control", "raw_v5_prior"])
        for index, sample_id in enumerate(ids):
            writer.writerow([
                sample_id,
                float(outputs["2d"][index]),
                float(outputs["real"][index]),
                float(outputs["shuffled"][index]),
                float(outputs["pocket"][index]),
                float(outputs["raw"][index]),
            ])

    print(json.dumps({
        "status": "ok",
        "model": "GeoFuse-CPI V2.0 native-structure expert",
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "n": n,
        "geometry_used": True,
        "input_contract": "precomputed legacy 423D v5 anchor + 704D native ligand-pocket descriptors",
        "architecture_note": "Separate frozen native residual expert; not jointly trained with the 2D V2.0 trunk.",
        "output": str(args.out),
        "notice": "Predictions are model estimates. This route does not generate poses or derive geometry from raw PDB/SDF files.",
    }, indent=2))


if __name__ == "__main__":
    main()
