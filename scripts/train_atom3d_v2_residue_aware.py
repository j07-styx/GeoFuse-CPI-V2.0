#!/usr/bin/env python3
"""Exploratory ATOM3D V2.0 residue-aware/physical 3D ablations.

Uses only the official 3,507-complex TRAIN split and frozen V2.0 ESM2/Morgan
representation. Nested 4-fold heuristic protein-cluster CV; official VAL/TEST
are not opened. All candidate weights/predictions stay local under the
ATOM3D CC BY-NC-ND restrictions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "inference"))
sys.path.insert(0, str(ROOT / "src"))
from geofuse_portable_kinase import Unified, unit
from geofuse_cpi.model import ExplicitPairResidual


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def group_weights(groups):
    _, inverse, counts = np.unique(groups, return_inverse=True, return_counts=True)
    weights = 1.0 / counts[inverse].astype(np.float64)
    return (weights / weights.mean()).astype(np.float32)


def fit_ridge(x, y, groups, alpha=1.0, balanced=False):
    mean, std = float(np.mean(y)), float(np.std(y))
    weights = group_weights(groups) if balanced else None
    model = Ridge(alpha=alpha, fit_intercept=True)
    model.fit(x, (y - mean) / std, sample_weight=weights)
    return model, mean, std


def predict_ridge(model, mean, std, x):
    return model.predict(x).astype(np.float64) * std + mean


def fit_residual(features, residual, y_scale, groups, seed, epochs=16):
    x = np.asarray(features, dtype=np.float32)
    y = np.asarray(residual, dtype=np.float32) / float(y_scale)
    mean = x.mean(0).astype(np.float32)
    std = x.std(0).astype(np.float32)
    std[std < 1e-6] = 1.0
    x = (x - mean) / std
    weights = group_weights(groups)

    torch.manual_seed(seed)
    net = ExplicitPairResidual(x.shape[1]).cpu()
    optimizer = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=2e-3)
    xt, yt, wt = torch.from_numpy(x), torch.from_numpy(y), torch.from_numpy(weights)
    rng = np.random.default_rng(seed)
    net.train()
    for epoch in range(epochs):
        order = rng.permutation(len(x))
        for start in range(0, len(order), 32):
            index = torch.from_numpy(order[start:start + 32].astype(np.int64))
            optimizer.zero_grad(set_to_none=True)
            error = net(xt[index]) - yt[index]
            loss = torch.mean(wt[index] * error.square())
            loss.backward()
            optimizer.step()
    net.eval()
    return {
        "state_dict": {k: v.detach().cpu().clone() for k, v in net.state_dict().items()},
        "mean": mean,
        "std": std,
        "scale": float(y_scale),
        "feature_dim": int(x.shape[1]),
    }


def residual_predict(expert, features):
    net = ExplicitPairResidual(expert["feature_dim"]).cpu()
    net.load_state_dict(expert["state_dict"])
    net.eval()
    x = (np.asarray(features, dtype=np.float32) - expert["mean"]) / expert["std"]
    with torch.inference_mode():
        return net(torch.from_numpy(x)).numpy().astype(np.float64) * expert["scale"]


def metrics(y, pred, group):
    group_rmse = [
        float(np.sqrt(np.mean((pred[group == g] - y[group == g]) ** 2)))
        for g in np.unique(group)
    ]
    return {
        "pooled_rmse": float(np.sqrt(mean_squared_error(y, pred))),
        "pooled_mae": float(mean_absolute_error(y, pred)),
        "pearson": float(np.corrcoef(y, pred)[0, 1]),
        "heuristic_protein_group_macro_rmse": float(np.mean(group_rmse)),
    }


def bootstrap_delta(y, a, b, group, seed, n_boot=5000):
    gids = np.unique(group)
    delta = []
    for g in gids:
        mask = group == g
        delta.append(
            np.sqrt(np.mean((a[mask] - y[mask]) ** 2))
            - np.sqrt(np.mean((b[mask] - y[mask]) ** 2))
        )
    delta = np.asarray(delta)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(gids), size=(n_boot, len(gids)))
    boot = delta[draws].mean(1)
    return {
        "mean_group_macro_delta_rmse": float(delta.mean()),
        "ci95": [float(v) for v in np.quantile(boot, [0.025, 0.975])],
        "groups": int(len(gids)),
        "replicates": int(n_boot),
        "negative_favors_first_arm": True,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-checkpoint", required=True, type=Path)
    ap.add_argument("--prepared", required=True, type=Path)
    ap.add_argument("--esm-embeddings", required=True, type=Path)
    ap.add_argument("--residue-features", required=True, type=Path)
    ap.add_argument("--physical-features", required=True, type=Path)
    ap.add_argument("--summary-out", required=True, type=Path)
    ap.add_argument("--predictions-out", required=True, type=Path)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--epochs", type=int, default=16)
    ap.add_argument("--bootstrap", type=int, default=5000)
    args = ap.parse_args()
    inputs = (args.base_checkpoint, args.prepared, args.esm_embeddings,
              args.residue_features, args.physical_features)
    if any(not p.is_file() for p in inputs):
        ap.error("A required local ATOM3D/V2.0 input is missing")
    if args.summary_out.exists() or args.predictions_out.exists():
        ap.error("Refusing to overwrite an existing result")
    torch.set_num_threads(4)

    prepared = np.load(args.prepared, allow_pickle=False)
    esm_cache = np.load(args.esm_embeddings, allow_pickle=False)
    residue = np.load(args.residue_features, allow_pickle=False)
    physical = np.load(args.physical_features, allow_pickle=False)
    y = prepared["y"].astype(np.float64)
    fold = prepared["fold_id"].astype(np.int8)
    group = prepared["protein_group_id"].astype(np.int32)
    if not np.array_equal(residue["y"], y):
        raise ValueError("Residue-aware feature labels/order do not match official TRAIN")
    if not np.array_equal(residue["original_train_indices"], prepared["source_index"]):
        raise ValueError("Residue-aware source indices do not match official TRAIN")
    if physical["real"].shape != (3507, 55):
        raise ValueError("Physical descriptor source order/count mismatch")
    if not all(np.isfinite(physical[k]).all() for k in ("real", "shuffled", "pocket")):
        raise ValueError("Nonfinite physical descriptor values")

    sample_index = esm_cache["sample_sequence_index"].astype(np.int32)
    protein = esm_cache["unique_embeddings"][sample_index].astype(np.float32)
    if not np.array_equal(
        esm_cache["unique_sequences"].astype(str)[sample_index],
        prepared["sequences"].astype(str),
    ):
        raise ValueError("ESM sequence mapping mismatch")
    compound = prepared["morgan"].astype(np.float32)
    base_ck = torch.load(args.base_checkpoint, map_location="cpu", weights_only=True)
    net = Unified().eval()
    missing, unexpected = net.load_state_dict(base_ck["state_dict"], strict=False)
    if unexpected or set(missing) - {"heads.atom3d.weight", "heads.atom3d.bias"}:
        raise ValueError(f"Base V2.0 model mismatch: {missing} / {unexpected}")
    p = unit(protein) * np.sqrt(1280)
    g = np.zeros((3507, 128), dtype=np.float32)
    latent = []
    with torch.inference_mode():
        for start in range(0, 3507, 256):
            stop = min(start + 256, 3507)
            latent.append(net.representation(
                torch.from_numpy(p[start:stop].astype(np.float32)),
                torch.from_numpy(compound[start:stop]),
                torch.from_numpy(g[start:stop]),
            ).numpy())
    latent = np.concatenate(latent).astype(np.float32)

    feature_sets = {
        "residue_1210": {k: residue[k].astype(np.float32) for k in ("real", "shuffled", "pocket")},
        "physical_55": {k: physical[k].astype(np.float32) for k in ("real", "shuffled", "pocket")},
        "combined_1265": {
            k: np.concatenate([residue[k], physical[k]], axis=1).astype(np.float32)
            for k in ("real", "shuffled", "pocket")
        },
    }
    oof_unweighted = np.full(3507, np.nan)
    oof_balanced = np.full(3507, np.nan)
    oof_3d = {source: {mode: np.full(3507, np.nan) for mode in ("real", "shuffled", "pocket")}
              for source in feature_sets}
    fold_records = []
    for outer in range(4):
        train = np.flatnonzero(fold != outer)
        test = np.flatnonzero(fold == outer)
        # Both 2D heads are evaluated on exactly the same outer fold.
        for balanced, target in ((False, oof_unweighted), (True, oof_balanced)):
            head, mu, sd = fit_ridge(
                latent[train], y[train], group[train], args.alpha, balanced)
            target[test] = predict_ridge(head, mu, sd, latent[test])

        # Nested OOF baseline predictions for residual training inside outer TRAIN.
        inner_unweighted = np.full(3507, np.nan)
        inner_balanced = np.full(3507, np.nan)
        for inner in range(4):
            if inner == outer:
                continue
            valid = np.flatnonzero((fold == inner) & (fold != outer))
            fit = np.flatnonzero((fold != inner) & (fold != outer))
            for balanced, target in ((False, inner_unweighted), (True, inner_balanced)):
                head, mu, sd = fit_ridge(
                    latent[fit], y[fit], group[fit], args.alpha, balanced)
                target[valid] = predict_ridge(head, mu, sd, latent[valid])
        if not np.isfinite(inner_balanced[train]).all():
            raise AssertionError("Nested group-balanced 2D OOF predictions incomplete")
        residual_train = y[train] - inner_balanced[train]
        scale = float(np.std(y[train]))
        record = {"outer_fold": outer, "n_train": int(len(train)), "n_test": int(len(test))}
        for source, arms in feature_sets.items():
            record[source] = {}
            for arm_i, mode in enumerate(("real", "shuffled", "pocket")):
                expert = fit_residual(
                    arms[mode][train], residual_train, scale, group[train],
                    seed=20261201 + outer, epochs=args.epochs)
                correction = residual_predict(expert, arms[mode][test])
                oof_3d[source][mode][test] = oof_balanced[test] + correction
                record[source][mode] = float(
                    np.sqrt(np.mean((oof_3d[source][mode][test] - y[test]) ** 2))
                )
        fold_records.append(record)
        print(f"Residue-aware ATOM3D outer fold {outer+1}/4 complete", flush=True)

    result_metrics = {
        "2d_unweighted": metrics(y, oof_unweighted, group),
        "2d_group_balanced": metrics(y, oof_balanced, group),
    }
    comparisons = {}
    for source in feature_sets:
        for mode in ("real", "shuffled", "pocket"):
            result_metrics[f"{source}_{mode}"] = metrics(y, oof_3d[source][mode], group)
        comparisons[source] = {
            "native_vs_group_balanced_2d": bootstrap_delta(
                y, oof_3d[source]["real"], oof_balanced, group, 20261211, args.bootstrap),
            "native_vs_shuffled": bootstrap_delta(
                y, oof_3d[source]["real"], oof_3d[source]["shuffled"], group,
                20261212, args.bootstrap),
            "native_vs_pocket": bootstrap_delta(
                y, oof_3d[source]["real"], oof_3d[source]["pocket"], group,
                20261213, args.bootstrap),
        }

    summary = {
        "status": "ATOM3D_V2_0_RESIDUE_AWARE_EXPLORATORY_CV_COMPLETE",
        "n_official_train": 3507,
        "official_validation_opened": False,
        "official_test_opened_or_rescored": False,
        "grouping": "1,004 observed CA-chain 3-mer Jaccard >=0.28 components; not aligned identity",
        "target": "pK from mixed Ki/Kd source records",
        "cv_metrics": result_metrics,
        "paired_group_bootstrap": comparisons,
        "outer_folds": fold_records,
        "license": "CC BY-NC-ND; local noncommercial research only; no redistribution of derivatives or weights without permission",
        "inputs_sha256": {
            "base_checkpoint": sha256(args.base_checkpoint),
            "prepared": sha256(args.prepared),
            "esm_embeddings": sha256(args.esm_embeddings),
            "residue_features": sha256(args.residue_features),
            "physical_features": sha256(args.physical_features),
        },
        "interpretation_limits": [
            "Earlier ATOM3D TRAIN labels informed V13-V36 exploratory work; these CV results are retrospective development evidence.",
            "The protein grouping is a 3-mer Jaccard proxy, not a true aligned sequence-identity split.",
            "The 55D physical descriptor is Vina-like and is not a calibrated docking energy or full AutoDock Vina.",
            "The pocket-only feature arm is selected using the experimental <=8 A contact graph and is a control, not a deployable pose-free feature source."
        ],
    }
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    np.savez_compressed(
        args.predictions_out, y=y.astype(np.float32), fold_id=fold, group_id=group,
        oof_2d_unweighted=oof_unweighted.astype(np.float32),
        oof_2d_group_balanced=oof_balanced.astype(np.float32),
        **{f"oof_{s}_{m}": p.astype(np.float32)
           for s, arms in oof_3d.items() for m, p in arms.items()},
    )
    print(json.dumps({"status": summary["status"], "cv_metrics": result_metrics,
                      "paired_group_bootstrap": comparisons}, indent=2), flush=True)


if __name__ == "__main__":
    main()
