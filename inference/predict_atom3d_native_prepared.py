#!/usr/bin/env python3
"""Prepared-complex inference for the ATOM3D-integrated GeoFuse-CPI V2.0 head.

This consumes a 2D sequence/SMILES input plus a precomputed 704D native
ligand-pocket descriptor. It does not dock ligands or generate poses.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "inference"))
sys.path.insert(0, str(ROOT / "src"))
from geofuse_portable_kinase import Predictor
from geofuse_cpi.model import ExplicitPairResidual


def predict_native(predictor, smiles, protein_embeddings, native_geometry,
                   allow_experimental=False):
    expert = predictor.ck.get("atom3d_native_expert")
    if not expert:
        raise ValueError("Checkpoint has no ATOM3D native-geometry expert")
    if not expert.get("deployment_eligible", False) and not allow_experimental:
        raise RuntimeError(
            "ATOM3D native expert is experimental: group-bootstrap evidence "
            "did not pass the prespecified deployment gate. Pass "
            "--allow-experimental only for research analysis."
        )
    geometry = np.asarray(native_geometry, dtype=np.float32)
    if geometry.shape != (len(smiles), 704) or not np.isfinite(geometry).all():
        raise ValueError("Expected one finite 704D native-geometry vector per pair")
    base = predictor.predict("atom3d", smiles, protein_embeddings)
    model = ExplicitPairResidual(704).cpu()
    model.load_state_dict(expert["state_dict"])
    model.eval()
    mean = expert["mean"]
    std = expert["std"]
    if isinstance(mean, torch.Tensor):
        mean = mean.detach().cpu().numpy()
    if isinstance(std, torch.Tensor):
        std = std.detach().cpu().numpy()
    standardized = (geometry - np.asarray(mean)) / np.asarray(std)
    with torch.inference_mode():
        correction = model(torch.from_numpy(standardized)).numpy()
    correction = correction * float(expert["scale"])
    return base, base + correction, correction


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--anchors", required=True, type=Path)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--esm-embeddings", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--allow-experimental", action="store_true")
    args = parser.parse_args()
    if args.out.exists():
        parser.error("Refusing to overwrite output")
    prepared = np.load(args.prepared, allow_pickle=False)
    esm_cache = np.load(args.esm_embeddings, allow_pickle=False)
    smiles = prepared["smiles"].astype(str)
    sample_index = esm_cache["sample_sequence_index"].astype(np.int32)
    embeddings = esm_cache["unique_embeddings"][sample_index].astype(np.float32)
    predictor = Predictor(args.checkpoint, args.anchors)
    base, native, correction = predict_native(
        predictor, smiles, embeddings, prepared["real"],
        allow_experimental=args.allow_experimental,
    )
    result = pd.DataFrame({
        "source_index": prepared["source_index"],
        "complex_id": prepared["ids"].astype(str),
        "predicted_pK_2D": base,
        "native_geometry_correction": correction,
        "predicted_pK_2D_plus_native_3D": native,
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)
    print(f"Saved {len(result)} prepared-complex predictions to {args.out}")


if __name__ == "__main__":
    main()
