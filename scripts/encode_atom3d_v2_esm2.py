#!/usr/bin/env python3
"""Compute frozen ESM-2 650M embeddings for ATOM3D official TRAIN only.

Outputs are local licensed-data derivatives and must not be redistributed.
Uses the exact V2.0 ESM layer-33 mean-residue embedding convention.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import torch
import esm


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_model(checkpoint: Path, device: torch.device):
    expected = "ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c"
    actual = sha256(checkpoint)
    if actual != expected:
        raise ValueError(f"Unexpected ESM-2 650M checkpoint SHA-256: {actual}")
    data = torch.load(checkpoint, map_location="cpu", weights_only=False, mmap=True)
    with torch.device("meta"):
        model, alphabet, state = esm.pretrained._load_model_and_alphabet_core_v2(data)
    report = model.load_state_dict(state, strict=False, assign=True)
    allowed_missing = {"contact_head.regression.weight", "contact_head.regression.bias"}
    if report.unexpected_keys or set(report.missing_keys) - allowed_missing:
        raise ValueError(f"Unexpected ESM state mismatch: {report}")
    if any(value.is_meta for _, value in model.named_buffers()):
        raise ValueError("ESM checkpoint left an unexpected meta buffer")
    if any(value.is_meta for name, value in model.named_parameters()
           if not name.startswith("contact_head.")):
        raise ValueError("ESM checkpoint left an unexpected meta parameter")
    # Contact regression is not used for layer-33 embeddings.
    model.contact_head = None
    model = model.eval().to(device)
    del state, data
    return model, alphabet.get_batch_converter()


def encode_one(model, converter, sequence: str, device: torch.device,
               window: int = 1000, stride: int = 872) -> np.ndarray:
    if not sequence or len(sequence) > 2048:
        raise ValueError(f"Unsupported sequence length: {len(sequence)}")
    accumulator = torch.zeros((len(sequence), 1280), dtype=torch.float32, device="cpu")
    counts = torch.zeros((len(sequence), 1), dtype=torch.float32, device="cpu")
    with torch.inference_mode():
        for start in range(0, len(sequence), stride):
            end = min(start + window, len(sequence))
            _, _, tokens = converter([("query", sequence[start:end])])
            result = model(tokens.to(device), repr_layers=[33], return_contacts=False)
            residue = result["representations"][33][0, 1:1 + end - start].float().cpu()
            accumulator[start:end] += residue
            counts[start:end] += 1.0
            del result, tokens, residue
    embedding = (accumulator / counts).half().float().mean(0).numpy()
    if embedding.shape != (1280,) or not np.isfinite(embedding).all():
        raise ValueError("Invalid ESM2 embedding")
    return embedding.astype(np.float32)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path, help="Final NPZ path")
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--max-sequences", type=int, default=0,
                        help="Smoke-test cap; zero means all unique sequences")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--window", type=int, default=1000)
    parser.add_argument("--stride", type=int, default=872)
    args = parser.parse_args()
    if not args.prepared.is_file() or not args.checkpoint.is_file():
        parser.error("Prepared TRAIN features or ESM checkpoint not found")
    if args.window < 128 or args.stride < 64 or args.stride > args.window:
        parser.error("Invalid ESM chunk window/stride")
    data = np.load(args.prepared, allow_pickle=False)
    if data["y"].shape != (3507,) or data["sequences"].shape != (3507,):
        raise ValueError("Unexpected prepared ATOM3D source archive")
    sequences = data["sequences"].astype(str).tolist()
    unique, sample_index = [], np.empty(len(sequences), dtype=np.int32)
    lookup = {}
    for i, sequence in enumerate(sequences):
        if sequence not in lookup:
            lookup[sequence] = len(unique)
            unique.append(sequence)
        sample_index[i] = lookup[sequence]
    if len(unique) != 2105:
        raise ValueError(f"Unexpected unique protein sequence count: {len(unique)}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    partial = args.out.with_suffix(".partial.npy")
    done_path = args.out.with_suffix(".done.npy")
    if args.resume:
        if not partial.exists() or not done_path.exists():
            raise FileNotFoundError("Resume requested but partial embedding files are missing")
        embeddings = np.load(partial, mmap_mode="r+")
        done = np.load(done_path)
        if embeddings.shape != (len(unique), 1280) or done.shape != (len(unique),):
            raise ValueError("Partial ESM cache shape mismatch")
    else:
        if args.out.exists() or partial.exists() or done_path.exists():
            raise FileExistsError("Output exists; use --resume or choose a new path")
        embeddings = np.lib.format.open_memmap(
            partial, mode="w+", dtype=np.float32, shape=(len(unique), 1280))
        done = np.zeros(len(unique), dtype=bool)
        np.save(done_path, done)

    device = torch.device(args.device)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable in this Python environment")
    torch.set_num_threads(2)
    if device.type == "cuda":
        torch.cuda.set_per_process_memory_fraction(0.82, device)
        torch.backends.cuda.matmul.allow_tf32 = False
    model, converter = load_model(args.checkpoint, device)
    limit = min(len(unique), args.max_sequences) if args.max_sequences else len(unique)
    for i in range(limit):
        if done[i]:
            continue
        try:
            embeddings[i] = encode_one(model, converter, unique[i], device,
                                       window=args.window, stride=args.stride)
        except torch.cuda.OutOfMemoryError:
            if device.type == "cuda":
                torch.cuda.empty_cache()
            raise
        done[i] = True
        if (i + 1) % 25 == 0 or i + 1 == limit:
            embeddings.flush()
            np.save(done_path, done)
            print(f"ESM2 ATOM3D unique sequences {i+1}/{limit}; completed {int(done.sum())}/{len(unique)}",
                  flush=True)
    embeddings.flush()
    np.save(done_path, done)
    if not bool(done.all()):
        print(f"PARTIAL: {int(done.sum())}/{len(unique)} embeddings complete. Re-run with --resume.",
              flush=True)
        return
    np.savez_compressed(args.out, unique_sequences=np.asarray(unique),
                        unique_embeddings=np.asarray(embeddings),
                        sample_sequence_index=sample_index)
    report = {
        "status": "ATOM3D_V2_0_ESM2_650M_TRAIN_ONLY_EMBEDDINGS",
        "n_complexes": 3507,
        "n_unique_sequences": len(unique),
        "embedding_shape": [len(unique), 1280],
        "esm_layer": 33,
        "esm_checkpoint_sha256": sha256(args.checkpoint),
        "prepared_train_sha256": sha256(args.prepared),
        "output_sha256": sha256(args.out),
        "official_val_or_test_opened": False,
        "official_test_rescored": False,
        "license": "CC BY-NC-ND; local noncommercial research only; do not redistribute dataset derivatives or trained weights without permission",
    }
    args.out.with_suffix(".audit.json").write_text(json.dumps(report, indent=2) + "\n",
                                                   encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
