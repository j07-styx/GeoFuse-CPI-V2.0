#!/usr/bin/env python3
"""Raw SMILES + FASTA inference for the frozen GeoFuse-CPI V2.0 2D predictor."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from geofuse_portable_kinase import Predictor
from geofuse_raw_kinase import Encoder, predict_raw

CHECKPOINT = ROOT / "models" / "v2.0" / "GEOFUSE_CPI_V2.pt"
ANCHORS = ROOT / "models" / "v2.0" / "TRAIN_ANCHORS_NO_Y.npz"
ESM_WEIGHTS = ROOT / "models" / "esm2_t33_650M_UR50D.pt"
CHECKPOINT_SHA256 = "b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4"
ANCHOR_SHA256 = "42ddd4d3af1f00f15bf96985adaaa15b821fb67fe6a5f2fddf3fa796a4773769"
UNITS = {"davis": "pKd", "kiba": "KIBA score", "metz": "pKd"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", required=True, help="One valid compound SMILES.")
    parser.add_argument("--fasta", required=True, type=Path, help="One protein FASTA file.")
    parser.add_argument("--assay", required=True, choices=("davis", "kiba", "metz"))
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT)
    parser.add_argument("--anchors", type=Path, default=ANCHORS)
    parser.add_argument("--esm-checkpoint", type=Path, default=ESM_WEIGHTS)
    args = parser.parse_args()

    if not args.fasta.is_file():
        parser.error(f"FASTA file not found: {args.fasta}")
    if not args.checkpoint.is_file() or not args.anchors.is_file():
        parser.error("Model assets missing. Run scripts/prepare_inference_assets.py.")
    if sha256_file(args.checkpoint) != CHECKPOINT_SHA256:
        parser.error("Checkpoint SHA-256 does not match the frozen V2.0 release.")
    if sha256_file(args.anchors) != ANCHOR_SHA256:
        parser.error("Training-anchor SHA-256 does not match the frozen V2.0 release.")
    if not args.esm_checkpoint.is_file():
        parser.error("ESM-2 weights missing. Run scripts/download_esm2_weights.py.")

    sequence = args.fasta.read_text(encoding="utf-8")
    predictor = Predictor(args.checkpoint, args.anchors)
    try:
        encoder = Encoder(args.esm_checkpoint)
        score = predict_raw(predictor, encoder, args.assay, args.smiles, sequence)
    finally:
        predictor.a.close()

    result = {
        "status": "ok",
        "model": "GeoFuse-CPI V2.0",
        "model_version": "geofuse-cpi-v2.0-frozen",
        "checkpoint_sha256": CHECKPOINT_SHA256,
        "assay_head": args.assay,
        "prediction": score,
        "unit": UNITS[args.assay],
        "geometry_used": False,
        "notice": "Model estimate, not an experimental measurement. This raw-input route is the 2D sequence-and-chemistry predictor; native geometry requires a separately prepared 423D/704D input.",
    }
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
