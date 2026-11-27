#!/usr/bin/env python3
"""Fit a group-balanced ATOM3D head + 704D native expert in GeoFuse-CPI V2.0.

Uses only official ATOM3D TRAIN (3,507 complexes), nested 4-fold protein
cluster CV, and a frozen V2.0 ESM2/Morgan trunk. No official VAL/TEST access.
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
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def group_weights(group):
    _, inverse, counts = np.unique(group, return_inverse=True, return_counts=True)
    weights = 1.0 / counts[inverse].astype(np.float64)
    return (weights / weights.mean()).astype(np.float32)


def fit_head(x, y, group, alpha):
    mean, std = float(y.mean()), float(y.std())
    model = Ridge(alpha=alpha, fit_intercept=True)
    model.fit(x, (y - mean) / std, sample_weight=group_weights(group))
    return model, mean, std


def head_predict(model, mean, std, x):
    return model.predict(x).astype(np.float64) * std + mean


def fit_residual(x, residual, y_scale, group, seed, epochs):
    x = np.asarray(x, dtype=np.float32)
    target = np.asarray(residual, dtype=np.float32) / float(y_scale)
    mean, std = x.mean(0).astype(np.float32), x.std(0).astype(np.float32)
    std[std < 1e-6] = 1.0
    x = (x - mean) / std
    weights = group_weights(group)
    torch.manual_seed(seed)
    model = ExplicitPairResidual(704).cpu()
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=2e-3)
    xt, yt, wt = torch.from_numpy(x), torch.from_numpy(target), torch.from_numpy(weights)
    rng = np.random.default_rng(seed)
    model.train()
    for epoch in range(epochs):
        order = rng.permutation(len(x))
        for start in range(0, len(order), 32):
            index = torch.from_numpy(order[start:start + 32].astype(np.int64))
            optimizer.zero_grad(set_to_none=True)
            error = model(xt[index]) - yt[index]
            loss = torch.mean(wt[index] * error.square())
            loss.backward()
            optimizer.step()
    model.eval()
    return {
        "state_dict": {k: v.detach().cpu().clone() for k, v in model.state_dict().items()},
        "mean": mean, "std": std, "scale": float(y_scale), "feature_dim": 704,
    }


def residual_predict(expert, x):
    model = ExplicitPairResidual(704).cpu()
    model.load_state_dict(expert["state_dict"])
    model.eval()
    values = (np.asarray(x, dtype=np.float32) - expert["mean"]) / expert["std"]
    with torch.inference_mode():
        return model(torch.from_numpy(values)).numpy().astype(np.float64) * expert["scale"]


def metrics(y, p, group):
    group_rmse = [
        np.sqrt(np.mean((p[group == g] - y[group == g]) ** 2))
        for g in np.unique(group)
    ]
    return {
        "pooled_rmse": float(np.sqrt(mean_squared_error(y, p))),
        "pooled_mae": float(mean_absolute_error(y, p)),
        "pearson": float(np.corrcoef(y, p)[0, 1]),
        "heuristic_protein_group_macro_rmse": float(np.mean(group_rmse)),
    }


def bootstrap_delta(y, a, b, group, seed, n_boot):
    gids = np.unique(group)
    delta = []
    for gid in gids:
        mask = group == gid
        delta.append(np.sqrt(np.mean((a[mask] - y[mask]) ** 2))
                     - np.sqrt(np.mean((b[mask] - y[mask]) ** 2)))
    delta = np.asarray(delta)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(gids), size=(n_boot, len(gids)))
    boot = delta[draws].mean(1)
    return {
        "mean_group_macro_delta_rmse": float(delta.mean()),
        "ci95": [float(v) for v in np.quantile(boot, [0.025, 0.975])],
        "groups": int(len(gids)), "replicates": int(n_boot),
        "negative_favors_first_arm": True,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base-checkpoint", required=True, type=Path)
    ap.add_argument("--prepared", required=True, type=Path)
    ap.add_argument("--esm-embeddings", required=True, type=Path)
    ap.add_argument("--out-checkpoint", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    ap.add_argument("--out-predictions", required=True, type=Path)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--epochs", type=int, default=16)
    ap.add_argument("--bootstrap", type=int, default=5000)
    args = ap.parse_args()
    for path in (args.base_checkpoint, args.prepared, args.esm_embeddings):
        if not path.is_file():
            ap.error(f"Missing input: {path}")
    if any(p.exists() for p in (args.out_checkpoint, args.out_summary, args.out_predictions)):
        ap.error("Refusing to overwrite an existing group-balanced candidate")
    torch.set_num_threads(4)

    prepared = np.load(args.prepared, allow_pickle=False)
    esm = np.load(args.esm_embeddings, allow_pickle=False)
    y = prepared["y"].astype(np.float64)
    group = prepared["protein_group_id"].astype(np.int32)
    fold = prepared["fold_id"].astype(np.int8)
    geometry = {k: prepared[k].astype(np.float32) for k in ("real", "shuffled", "pocket")}
    sample_index = esm["sample_sequence_index"].astype(np.int32)
    protein = esm["unique_embeddings"][sample_index].astype(np.float32)
    if not np.array_equal(esm["unique_sequences"].astype(str)[sample_index],
                          prepared["sequences"].astype(str)):
        raise ValueError("ESM sequence mapping mismatch")
    compound = prepared["morgan"].astype(np.float32)

    base = torch.load(args.base_checkpoint, map_location="cpu", weights_only=True)
    model = Unified().eval()
    missing, unexpected = model.load_state_dict(base["state_dict"], strict=False)
    if unexpected or set(missing) - {"heads.atom3d.weight", "heads.atom3d.bias"}:
        raise ValueError(f"V2.0 checkpoint mismatch: {missing}/{unexpected}")
    p = unit(protein) * np.sqrt(1280)
    g = np.zeros((3507, 128), dtype=np.float32)
    latent = []
    with torch.inference_mode():
        for start in range(0, 3507, 256):
            stop = min(start + 256, 3507)
            latent.append(model.representation(
                torch.from_numpy(p[start:stop].astype(np.float32)),
                torch.from_numpy(compound[start:stop]),
                torch.from_numpy(g[start:stop]),
            ).numpy())
    latent = np.concatenate(latent).astype(np.float32)

    oof_2d = np.full(3507, np.nan)
    oof_3d = {mode: np.full(3507, np.nan) for mode in geometry}
    fold_records = []
    for outer in range(4):
        train, test = np.flatnonzero(fold != outer), np.flatnonzero(fold == outer)
        outer_head, mu, sd = fit_head(latent[train], y[train], group[train], args.alpha)
        base_test = head_predict(outer_head, mu, sd, latent[test])
        oof_2d[test] = base_test

        inner_oof = np.full(3507, np.nan)
        for inner in range(4):
            if inner == outer:
                continue
            valid = np.flatnonzero((fold == inner) & (fold != outer))
            fit = np.flatnonzero((fold != inner) & (fold != outer))
            inner_head, im, istd = fit_head(latent[fit], y[fit], group[fit], args.alpha)
            inner_oof[valid] = head_predict(inner_head, im, istd, latent[valid])
        if not np.isfinite(inner_oof[train]).all():
            raise AssertionError("Nested group-balanced OOF head incomplete")
        residual = y[train] - inner_oof[train]
        record = {"outer_fold": outer, "n_train": len(train), "n_test": len(test)}
        for mode in ("real", "shuffled", "pocket"):
            expert = fit_residual(
                geometry[mode][train], residual, float(np.std(y[train])),
                group[train], seed=20261220 + outer, epochs=args.epochs)
            oof_3d[mode][test] = base_test + residual_predict(expert, geometry[mode][test])
            record[mode] = float(np.sqrt(np.mean((oof_3d[mode][test] - y[test]) ** 2)))
        fold_records.append(record)
        print(f"Group-balanced V2.0 704D ATOM3D outer fold {outer+1}/4 complete", flush=True)

    cv = {"2d_group_balanced": metrics(y, oof_2d, group)}
    for mode in oof_3d:
        cv[mode] = metrics(y, oof_3d[mode], group)
    paired = {
        "native_vs_2d": bootstrap_delta(y, oof_3d["real"], oof_2d, group, 20261221, args.bootstrap),
        "native_vs_shuffled": bootstrap_delta(y, oof_3d["real"], oof_3d["shuffled"], group, 20261222, args.bootstrap),
        "native_vs_pocket": bootstrap_delta(y, oof_3d["real"], oof_3d["pocket"], group, 20261223, args.bootstrap),
    }
    eligible = (
        cv["real"]["pooled_rmse"] < cv["2d_group_balanced"]["pooled_rmse"]
        and cv["real"]["heuristic_protein_group_macro_rmse"]
            < cv["2d_group_balanced"]["heuristic_protein_group_macro_rmse"]
        and paired["native_vs_2d"]["ci95"][1] < 0
        and paired["native_vs_shuffled"]["ci95"][1] < 0
        and paired["native_vs_pocket"]["ci95"][1] < 0
    )

    # Final group-balanced 2D task head.
    final_mean, final_std = float(y.mean()), float(y.std())
    weights = group_weights(group)
    final_head = Ridge(alpha=args.alpha, fit_intercept=True)
    final_head.fit(latent, (y - final_mean) / final_std, sample_weight=weights)
    with torch.no_grad():
        model.heads.atom3d.weight.copy_(
            torch.from_numpy(final_head.coef_.astype(np.float32).reshape(1, 48)))
        model.heads.atom3d.bias.copy_(torch.tensor([float(final_head.intercept_)]))
    base["state_dict"] = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    base["assay_normalization"] = dict(base["assay_normalization"])
    base["assay_normalization"]["atom3d"] = [final_mean, final_std]

    final_residual = y - oof_2d
    final_expert = fit_residual(
        geometry["real"], final_residual, float(np.std(y)), group,
        seed=20261227, epochs=args.epochs)
    base["atom3d_native_expert"] = {
        "state_dict": final_expert["state_dict"],
        "mean": torch.from_numpy(final_expert["mean"].copy()),
        "std": torch.from_numpy(final_expert["std"].copy()),
        "scale": final_expert["scale"],
        "feature_dim": 704,
        "training_n": 3507,
        "target": "residual to nested out-of-fold group-balanced V2.0 ATOM3D 2D head",
        "geometry": "704D experimentally bound native ligand-pocket descriptor",
        "status": "ELIGIBLE" if eligible else "EXPERIMENTAL_NOT_STATISTICALLY_ESTABLISHED",
        "deployment_eligible": bool(eligible),
    }
    base["atom3d_training_metadata"] = dict(base.get("atom3d_training_metadata", {}))
    base["atom3d_training_metadata"].update({
        "2d_head_training": "Group-balanced Ridge alpha=1 over frozen V2.0 fused 48D representation",
        "native_geometry_training": "Group-balanced 704D residual expert trained on nested OOF residuals",
        "cv_metrics_group_balanced": cv,
        "paired_group_bootstrap_group_balanced": paired,
        "native_geometry_deployment_eligible": bool(eligible),
    })

    args.out_checkpoint.parent.mkdir(parents=True, exist_ok=True)
    torch.save(base, args.out_checkpoint)
    np.savez_compressed(args.out_predictions, y=y.astype(np.float32), fold_id=fold,
                        group_id=group, oof_2d=oof_2d.astype(np.float32),
                        **{f"oof_{k}": v.astype(np.float32) for k, v in oof_3d.items()})
    summary = {
        "status": "ATOM3D_V2_0_GROUP_BALANCED_2D_AND_704D_NATIVE_EXPERT_TRAINED",
        "n_official_train": 3507,
        "official_validation_opened": False,
        "official_test_opened_or_rescored": False,
        "cv_metrics": cv,
        "paired_group_bootstrap": paired,
        "native_geometry_deployment_eligible": bool(eligible),
        "checkpoint_sha256": sha256(args.out_checkpoint),
        "predictions_sha256": sha256(args.out_predictions),
        "existing_davis_kiba_metz_gensparc_weights_unchanged": True,
        "license": "CC BY-NC-ND; local noncommercial research only; no redistribution without permission",
        "outer_folds": fold_records,
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
