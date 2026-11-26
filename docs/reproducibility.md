# Reproducibility and verification

## What is public in this repository

This repository is a **documentation and aggregate-evidence record**. It includes methodology and architecture documentation, aggregate-only benchmark and native-geometry JSON records, source/checkpoint identity hashes, explicit validation gates and limitations, and a metadata consistency checker.

It does **not** distribute model checkpoint binaries, raw affinity datasets, protected labels, per-example prediction arrays, raw PDB/SDF structures, ESM-2 weights or private source-training records.

## Frozen artifact identity

- Model checkpoint SHA-256: `b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4`
- Prediction archive SHA-256: `b8f935e2115ead2b83e247477118adedd38559acf368b80eabafc160240f0555`
- The original local V2.0 release manifest pins files in the local inference bundle. The bundle itself is not redistributed here.

## Public metadata check

From the repository root:

```bash
python scripts/verify_public_records.py
```

This checks the internal consistency of the public JSON records only. It does not run model inference or verify a checkpoint that is not included.

## Original local release verification

On the original research machine, the following checks were run successfully using the Python 3.11 environment and the local frozen V2.0 release bundle:

```bash
cd research_v2
python verify_release.py
python verify_release.py --smoke
```

Both completed with: `PASS: one checkpoint, label-free bundle and result ledger verified.`

The smoke check verifies bundle hashes and finite example predictions. It is **not** a new accuracy evaluation, a retraining run, or proof of clean-room installation on a different computer.

## Retrospective native replay

The original release's replay instructions use `python reproduce_biolip3.py` and `python reproduce.py` with the original pinned replay dependencies and authorized local source snapshots. Those scripts recompute the fourfold model outputs and compare them with frozen references before reading source affinities. They do not retrain the V2.0 2D backbone, regenerate every raw ESM-2/native feature from first principles, retrain PerceiverCPI, or reproduce the official fivefold DAVIS/KIBA/METZ protocols.

The complete replay bundles and scripts remain in the original local research checkout and are not included in this documentation-only repository. A reader cannot reproduce all model predictions from a fresh clone of this repository alone.

## Reproduction levels

1. **Public-record verification:** available from this repository.
2. **Frozen-bundle integrity and smoke inference:** passed on the original machine; requires the local release bundle and dependencies.
3. **Retrospective BioLiP3/OpenBind replay:** documented in the original local release; requires permitted source data snapshots and replay bundles.
4. **Full independent retraining / clean-room reproduction:** not established by the current public release.

Do not describe level 1 or 2 as full independent reproduction.
