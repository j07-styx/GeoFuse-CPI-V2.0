#!/usr/bin/env python3
"""Verify frozen V2.0 assets, prepared-input prediction, and optional raw-input parity."""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "inference"))

from geofuse_portable_kinase import Predictor
from geofuse_raw_kinase import Encoder, predict_raw

CHECKPOINT_SHA = "b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4"
ANCHOR_SHA = "42ddd4d3af1f00f15bf96985adaaa15b821fb67fe6a5f2fddf3fa796a4773769"
ESM_SHA = "ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c"
EXPECTED_ASPIRIN_AAK1 = 5.476265008935864


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", action="store_true", help="Run the full ESM-2 raw FASTA + SMILES smoke test.")
    parser.add_argument("--replay", action="store_true", help="Run exact OpenBind and BioLiP3 prepared-prediction parity.")
    parser.add_argument("--esm-checkpoint", type=Path, default=ROOT / "models" / "esm2_t33_650M_UR50D.pt")
    args = parser.parse_args()

    model_dir = ROOT / "models" / "v2.0"
    checkpoint = model_dir / "GEOFUSE_CPI_V2.pt"
    anchors = model_dir / "TRAIN_ANCHORS_NO_Y.npz"
    if not checkpoint.is_file() or not anchors.is_file():
        raise FileNotFoundError("Run scripts/prepare_inference_assets.py first.")
    if sha256(checkpoint) != CHECKPOINT_SHA or sha256(anchors) != ANCHOR_SHA:
        raise RuntimeError("Frozen V2.0 checkpoint or anchor hash mismatch.")

    predictor = Predictor(checkpoint, anchors)
    try:
        embedding = np.random.default_rng(2931).standard_normal((1, 1280)).astype(np.float32)
        smoke = predictor.predict("davis", ["CC(=O)Oc1ccccc1C(=O)O"], embedding)
        if smoke.shape != (1,) or not np.isfinite(smoke).all():
            raise RuntimeError("Prepared-embedding smoke prediction failed.")
        print(f"PASS: prepared-input DAVIS smoke prediction is finite ({smoke[0]:.6f}).")
    finally:
        predictor.a.close()

    if args.raw:
        if not args.esm_checkpoint.is_file():
            raise FileNotFoundError("ESM weights missing; run scripts/download_esm2_weights.py.")
        if sha256(args.esm_checkpoint) != ESM_SHA:
            raise RuntimeError("ESM-2 checkpoint hash mismatch.")
        sequence = (ROOT / "examples" / "AAK1.fasta").read_text(encoding="utf-8")
        predictor = Predictor(checkpoint, anchors)
        try:
            score = predict_raw(predictor, Encoder(args.esm_checkpoint), "davis",
                                "CC(=O)Oc1ccccc1C(=O)O", sequence)
        finally:
            predictor.a.close()
        if abs(score - EXPECTED_ASPIRIN_AAK1) > 1e-3:
            raise RuntimeError(f"Raw-input smoke parity failed: {score:.9f}")
        print(f"PASS: raw SMILES + FASTA parity ({score:.6f} pKd; tolerance 0.001).")

    if args.replay:
        commands = [
            [sys.executable, str(ROOT / "reproduce.py"), "--inference-only"],
            [sys.executable, str(ROOT / "reproduce_biolip3.py"), "--inference-only"],
        ]
        for command in commands:
            subprocess.run(command, cwd=ROOT, check=True)
        print("PASS: OpenBind (660) and BioLiP3 (73) frozen native prediction parity.")

    print("PASS: frozen V2.0 inference verification completed.")


if __name__ == "__main__":
    main()
