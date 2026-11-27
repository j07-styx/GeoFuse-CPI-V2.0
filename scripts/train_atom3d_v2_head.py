#!/usr/bin/env python3
"""Train the ATOM3D task head and native-geometry residual inside GeoFuse-CPI V2.0.

Only the official 3,507-complex TRAIN LMDB is used. Evaluation is nested,
heuristic sequence-cluster-held-out CV within TRAIN. Official VAL/TEST are not
opened. The resulting checkpoint is a local research artifact under the
ATOM3D CC BY-NC-ND restrictions; do not redistribute without permission.
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
from geofuse_portable_kinase import Unified, fingerprints, unit
from geofuse_cpi.model import ExplicitPairResidual


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def metrics(y: np.ndarray, pred: np.ndarray, group: np.ndarray) -> dict:
    pearson = float(np.corrcoef(y, pred)[0, 1])
    groups = np.unique(group)
    group_rmse = [
        float(np.sqrt(np.mean(np.square(pred[group == g] - y[group == g]))))
        for g in groups
    ]
    return {
        "pooled_rmse": float(np.sqrt(mean_squared_error(y, pred))),
        "pooled_mae": float(mean_absolute_error(y, pred)),
        "pearson": pearson,
        "heuristic_protein_group_macro_rmse": float(np.mean(group_rmse)),
        "n": int(len(y)),
        "n_heuristic_protein_groups": int(len(groups)),
    }


def group_bootstrap_macro_delta(y, pred_a, pred_b, group, seed=20261127, n_boot=5000):
    unique = np.unique(group)
    delta = np.empty(len(unique), dtype=np.float64)
    for j, gid in enumerate(unique):
        mask = group == gid
        rmse_a = np.sqrt(np.mean(np.square(pred_a[mask] - y[mask])))
        rmse_b = np.sqrt(np.mean(np.square(pred_b[mask] - y[mask])))
        delta[j] = rmse_a - rmse_b
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(unique), size=(n_boot, len(unique)))
    boot = delta[draws].mean(axis=1)
    return {
        "mean_group_macro_delta_rmse": float(delta.mean()),
        "ci95": [float(x) for x in np.quantile(boot, [0.025, 0.975])],
        "bootstrap_groups": int(len(unique)),
        "bootstrap_replicates": int(n_boot),
        "interpretation": "Negative means the first prediction arm has lower group-macro RMSE.",
    }


def fit_ridge(x_train, y_train, alpha=1.0):
    mean = float(np.mean(y_train))
    std = float(np.std(y_train))
    if std < 1e-8:
        raise ValueError("Degenerate ATOM3D training target scale")
    model = Ridge(alpha=alpha, fit_intercept=True)
    model.fit(x_train, (y_train - mean) / std)
    return model, mean, std


def ridge_predict(model, mean, std, x):
    return model.predict(x).astype(np.float64) * std + mean


def fit_residual(features, residual, y_scale, seed, epochs=16, batch_size=32):
    x = np.asarray(features, dtype=np.float32)
    y = np.asarray(residual, dtype=np.float32) / float(y_scale)
    mean = x.mean(axis=0).astype(np.float32)
    std = x.std(axis=0).astype(np.float32)
    std[std < 1e-6] = 1.0
    x = (x - mean) / std

    torch.manual_seed(seed)
    net = ExplicitPairResidual(704).cpu()
    optimizer = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=2e-3)
    xt = torch.from_numpy(x)
    yt = torch.from_numpy(y)
    rng = np.random.default_rng(seed)
    net.train()
    for epoch in range(epochs):
        order = rng.permutation(len(x))
        for start in range(0, len(order), batch_size):
            index = torch.from_numpy(order[start:start + batch_size].astype(np.int64))
            optimizer.zero_grad(set_to_none=True)
            prediction = net(xt[index])
            loss = F.mse_loss(prediction, yt[index])
            loss.backward()
            optimizer.step()
    net.eval()
    return {
        "state_dict": {k: v.detach().cpu().clone() for k, v in net.state_dict().items()},
        "mean": mean,
        "std": std,
        "scale": float(y_scale),
    }


def residual_predict(expert, features):
    net = ExplicitPairResidual(704).cpu()
    net.load_state_dict(expert["state_dict"])
    net.eval()
    x = (np.asarray(features, dtype=np.float32) - expert["mean"]) / expert["std"]
    with torch.inference_mode():
        correction = net(torch.from_numpy(x)).numpy().astype(np.float64)
    return correction * float(expert["scale"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-checkpoint", required=True, type=Path)
    parser.add_argument("--anchors", required=True, type=Path)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--esm-embeddings", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--predictions-out", required=True, type=Path)
    parser.add_argument("--summary-out", required=True, type=Path)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--epochs", type=int, default=16)
    parser.add_argument("--bootstrap", type=int, default=5000)
    args = parser.parse_args()
    for path in (args.base_checkpoint, args.anchors, args.prepared, args.esm_embeddings):
        if not path.is_file():
            parser.error(f"Required input missing: {path}")
    if any(path.exists() for path in (args.out, args.predictions_out, args.summary_out)):
        parser.error("Refusing to overwrite an existing V2.0 integration artifact")
    torch.set_num_threads(4)
    torch.manual_seed(20261127)

    prepared = np.load(args.prepared, allow_pickle=False)
    esm_cache = np.load(args.esm_embeddings, allow_pickle=False)
    y = prepared["y"].astype(np.float64)
    smiles = prepared["smiles"].astype(str)
    compound = prepared["morgan"].astype(np.float32)
    geometry = {
        mode: prepared[mode].astype(np.float32)
        for mode in ("real", "shuffled", "pocket")
    }
    fold_id = prepared["fold_id"].astype(np.int8)
    group_id = prepared["protein_group_id"].astype(np.int32)
    sample_sequence_index = esm_cache["sample_sequence_index"].astype(np.int32)
    unique_sequences = esm_cache["unique_sequences"].astype(str)
    prepared_sequences = prepared["sequences"].astype(str)
    if not np.array_equal(unique_sequences[sample_sequence_index], prepared_sequences):
        raise ValueError("ESM cache sequence-to-sample mapping does not match prepared TRAIN data")
    protein = esm_cache["unique_embeddings"][sample_sequence_index].astype(np.float32)
    if (len(y) != 3507 or protein.shape != (3507, 1280)
            or compound.shape != (3507, 2048)
            or any(value.shape != (3507, 704) for value in geometry.values())):
        raise ValueError("Unexpected ATOM3D V2.0 input shapes")
    if not all(np.isfinite(value).all() for value in (y, protein, compound, *geometry.values())):
        raise ValueError("Nonfinite ATOM3D training features")
    if set(np.unique(fold_id)) != {0, 1, 2, 3}:
        raise ValueError("Expected four frozen protein-cluster folds")
    if len(np.unique(group_id)) != 1004:
        raise ValueError("Expected 1,004 heuristic protein groups")
    for gid in np.unique(group_id):
        if len(np.unique(fold_id[group_id == gid])) != 1:
            raise ValueError("Heuristic protein group crosses CV folds")

    # Exact feature parity with the raw V2.0 Morgan inference helper.
    recomputed_fp = fingerprints(smiles)
    if not np.array_equal(recomputed_fp, compound):
        mismatch = int(np.sum(np.any(recomputed_fp != compound, axis=1)))
        raise ValueError(f"Morgan inference/training parity failed on {mismatch} rows")

    base = torch.load(args.base_checkpoint, map_location="cpu", weights_only=True)
    net = Unified()
    missing, unexpected = net.load_state_dict(base["state_dict"], strict=False)
    if unexpected or set(missing) - {"heads.atom3d.weight", "heads.atom3d.bias"}:
        raise ValueError(f"Base V2.0 checkpoint mismatch: missing={missing}, unexpected={unexpected}")
    net.eval()
    protein_input = unit(protein) * np.sqrt(1280)
    zeros = np.zeros((len(y), 128), dtype=np.float32)
    latent_rows = []
    with torch.inference_mode():
        for start in range(0, len(y), 256):
            stop = min(start + 256, len(y))
            latent_rows.append(net.representation(
                torch.from_numpy(protein_input[start:stop].astype(np.float32)),
                torch.from_numpy(compound[start:stop]),
                torch.from_numpy(zeros[start:stop]),
            ).numpy())
    latent = np.concatenate(latent_rows).astype(np.float32)
    if latent.shape != (3507, 48) or not np.isfinite(latent).all():
        raise ValueError("V2.0 shared latent representation is invalid")

    # Nested group-fold evaluation: each outer test fold is excluded from both
    # the 2D head and all residual-expert fitting.
    oof_2d = np.full(len(y), np.nan, dtype=np.float64)
    oof_3d = {mode: np.full(len(y), np.nan, dtype=np.float64) for mode in geometry}
    fold_results = []
    for outer in range(4):
        train = np.flatnonzero(fold_id != outer)
        test = np.flatnonzero(fold_id == outer)
        outer_head, outer_mean, outer_std = fit_ridge(latent[train], y[train], args.alpha)
        base_test = ridge_predict(outer_head, outer_mean, outer_std, latent[test])
        oof_2d[test] = base_test

        inner_oof = np.full(len(y), np.nan, dtype=np.float64)
        for inner in range(4):
            if inner == outer:
                continue
            valid = np.flatnonzero((fold_id == inner) & (fold_id != outer))
            fit = np.flatnonzero((fold_id != inner) & (fold_id != outer))
            inner_head, inner_mean, inner_std = fit_ridge(latent[fit], y[fit], args.alpha)
            inner_oof[valid] = ridge_predict(inner_head, inner_mean, inner_std, latent[valid])
        if not np.isfinite(inner_oof[train]).all():
            raise AssertionError("Nested 2D OOF predictions are incomplete")

        residual_train = y[train] - inner_oof[train]
        y_scale = float(np.std(y[train]))
        fold_record = {
            "outer_fold": outer,
            "n_train": int(len(train)),
            "n_test": int(len(test)),
            "2d_rmse": float(np.sqrt(np.mean(np.square(base_test - y[test])))),
            "3d_arms": {},
        }
        for arm_index, mode in enumerate(("real", "shuffled", "pocket")):
            expert = fit_residual(
                geometry[mode][train], residual_train, y_scale,
                seed=20261127 + outer, epochs=args.epochs,
            )
            correction = residual_predict(expert, geometry[mode][test])
            oof_3d[mode][test] = base_test + correction
            fold_record["3d_arms"][mode] = {
                "rmse": float(np.sqrt(np.mean(np.square(oof_3d[mode][test] - y[test]))))
            }
        fold_results.append(fold_record)
        print(f"ATOM3D V2.0 nested CV outer fold {outer+1}/4 complete", flush=True)

    if not np.isfinite(oof_2d).all() or any(not np.isfinite(v).all() for v in oof_3d.values()):
        raise AssertionError("Nested OOF predictions are incomplete")

    cv_metrics = {"2d": metrics(y, oof_2d, group_id)}
    for mode, prediction in oof_3d.items():
        cv_metrics[mode] = metrics(y, prediction, group_id)
    comparisons = {
        "native_vs_2d": group_bootstrap_macro_delta(
            y, oof_3d["real"], oof_2d, group_id, n_boot=args.bootstrap),
        "native_vs_shuffled": group_bootstrap_macro_delta(
            y, oof_3d["real"], oof_3d["shuffled"], group_id,
            seed=20261128, n_boot=args.bootstrap),
        "native_vs_pocket": group_bootstrap_macro_delta(
            y, oof_3d["real"], oof_3d["pocket"], group_id,
            seed=20261129, n_boot=args.bootstrap),
    }
    native_eligible = (
        cv_metrics["real"]["pooled_rmse"] < cv_metrics["2d"]["pooled_rmse"]
        and cv_metrics["real"]["heuristic_protein_group_macro_rmse"]
            < cv_metrics["2d"]["heuristic_protein_group_macro_rmse"]
        and comparisons["native_vs_2d"]["ci95"][1] < 0
        and comparisons["native_vs_shuffled"]["ci95"][1] < 0
        and comparisons["native_vs_pocket"]["ci95"][1] < 0
    )

    # Final 2D task head: frozen shared V2.0 tower, new linear ATOM3D pK head.
    final_mean = float(y.mean())
    final_std = float(y.std())
    final_head = Ridge(alpha=args.alpha, fit_intercept=True)
    final_head.fit(latent, (y - final_mean) / final_std)
    with torch.no_grad():
        net.heads.atom3d.weight.copy_(
            torch.from_numpy(final_head.coef_.astype(np.float32).reshape(1, 48)))
        net.heads.atom3d.bias.copy_(
            torch.tensor([float(final_head.intercept_)], dtype=torch.float32))
    base["state_dict"] = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
    base["assay_normalization"] = dict(base["assay_normalization"])
    base["assay_normalization"]["atom3d"] = [final_mean, final_std]

    # Fit one final native residual expert on OOF 2D residual targets. It is
    # stored even when the prespecified eligibility gate fails, but inference
    # must treat it as experimental unless the gate passes.
    final_residual = y - oof_2d
    final_native = fit_residual(
        geometry["real"], final_residual, float(np.std(y)),
        seed=20261227, epochs=args.epochs,
    )
    base["atom3d_native_expert"] = {
        "state_dict": final_native["state_dict"],
        "mean": torch.from_numpy(final_native["mean"].copy()),
        "std": torch.from_numpy(final_native["std"].copy()),
        "scale": final_native["scale"],
        "feature_dim": 704,
        "training_n": 3507,
        "target": "residual to nested out-of-fold V2.0 ATOM3D 2D head",
        "geometry": "experimentally bound native ligand-pocket coordinates",
        "status": "ELIGIBLE" if native_eligible else "EXPERIMENTAL_NOT_STATISTICALLY_ESTABLISHED",
        "deployment_eligible": bool(native_eligible),
    }
    base["atom3d_training_metadata"] = {
        "dataset": "ATOM3D LBA 30% sequence-identity split, official TRAIN only",
        "n_train": 3507,
        "n_unique_protein_sequences": int(len(unique_sequences)),
        "n_heuristic_protein_groups": int(len(np.unique(group_id))),
        "fold_sizes": [int(np.sum(fold_id == i)) for i in range(4)],
        "grouping": "Observed CA-chain 3-mer Jaccard >=0.28 connected components; not aligned percent identity",
        "target": "pK = -log10(K), mixed Ki/Kd source records",
        "shared_v2_0_tower_frozen": True,
        "new_task_head": "Ridge-fitted linear head over the frozen 48D V2.0 fused representation",
        "native_geometry_residual": "704D V2.0-compatible pair-distance/bond-angle descriptor; nested OOF 2D residual target",
        "official_validation_opened": False,
        "official_test_opened_or_rescored": False,
        "historical_exposure": "ATOM3D TRAIN labels were used in earlier exploratory V13-V36 branches; this is retrospective development, not independent validation",
        "license": "CC BY-NC-ND; local noncommercial research only; do not redistribute trained weights or derivatives without permission",
        "base_checkpoint_sha256": sha256(args.base_checkpoint),
        "prepared_train_features_sha256": sha256(args.prepared),
        "esm_embeddings_sha256": sha256(args.esm_embeddings),
        "cv_metrics": cv_metrics,
        "paired_group_bootstrap": comparisons,
        "native_geometry_deployment_eligible": bool(native_eligible),
        "outer_fold_results": fold_results,
        "alpha": float(args.alpha),
        "residual_epochs": int(args.epochs),
        "bootstrap_replicates": int(args.bootstrap),
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    torch.save(base, args.out)
    prediction_archive = {
        "y": y.astype(np.float32),
        "fold_id": fold_id,
        "protein_group_id": group_id,
        "oof_2d": oof_2d.astype(np.float32),
        **{f"oof_{mode}": values.astype(np.float32) for mode, values in oof_3d.items()},
    }
    np.savez_compressed(args.predictions_out, **prediction_archive)
    summary = {
        "status": "ATOM3D_TRAIN_ONLY_V2_0_TASK_HEAD_AND_NATIVE_EXPERT_TRAINED",
        "model_version": "GeoFuse-CPI V2.0",
        "checkpoint": str(args.out),
        "checkpoint_sha256": sha256(args.out),
        "n_official_train": 3507,
        "official_validation_opened": False,
        "official_test_opened_or_rescored": False,
        "cv_metrics": cv_metrics,
        "paired_group_bootstrap": comparisons,
        "native_geometry_deployment_eligible": bool(native_eligible),
        "existing_v2_0_heads_unchanged": True,
        "native_expert_included_as_experimental_if_gate_fails": True,
        "license": "CC BY-NC-ND; local noncommercial research only; do not redistribute trained weights or derivatives without permission",
        "prediction_archive_sha256": sha256(args.predictions_out),
    }
    args.summary_out.parent.mkdir(parents=True, exist_ok=True)
    args.summary_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
