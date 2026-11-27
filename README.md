# GeoFuse-CPI V2.0

**Frozen research release for GeoFuse-CPI V2.0.** This repository records the single frozen hybrid model, its architecture, benchmark evidence, provenance, inference implementation and reproducibility boundaries. It is separate from the production website repository and legacy experimental branches.

> **Release status:** frozen inference release candidate. The published checkpoint can now be run from raw SMILES + protein FASTA for the 2D predictor. Native-structure inference is available for prepared 423D/704D features. Full independent retraining and clean-room reproduction of every paper table are not yet established.

## Frozen artifact identity

- Checkpoint SHA-256: b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4
- Training-anchor SHA-256: 42ddd4d3af1f00f15bf96985adaaa15b821fb67fe6a5f2fddf3fa796a4773769
- Inference asset archive SHA-256: 71148521ae43187fcce4183cd7d388dc79daa04ebd6ace98bfacc08685a37017
- ESM-2 650M weights are downloaded separately and checked against SHA-256 ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c.
- One frozen checkpoint is used for the V2.0 headline ledger. No V363/V499/v6/V48 results are mixed into the 2D headline results.

## Quick start: raw 2D inference

Use Python 3.11.9. The tested reference environment is Windows CPU.

    python -m pip install torch==2.14.0
    python -m pip install -r inference/requirements-inference.txt
    python scripts/prepare_inference_assets.py
    python scripts/download_esm2_weights.py
    python inference/predict.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --fasta examples/AAK1.fasta --assay davis

The command returns a model estimate in the selected assay's scale. DAVIS and METZ are reported as pKd; KIBA is reported on its own KIBA scale. This is not an experimental measurement.

To run the bundled software smoke test:

    python scripts/verify_frozen_inference.py --raw

## Native-structure inference

The frozen checkpoint also contains the separate native-structure residual expert. It accepts **prepared, label-free features**, not raw PDB/SDF structures:

- x2d: (N, 423) legacy V5 2D anchor context
- real, shuffled, pocket: each (N, 704) geometry descriptors
- optional unique ids

    python inference/predict_native_prepared.py --features path/to/prepared_noY.npz --out predictions.csv

The original public native expert is not jointly trained with the V2.0 48-dimensional 2D trunk; it uses the legacy 423D V5 anchor. A separate local ATOM3D extension now adds a fifth task head on the frozen V2.0 shared representation and a distinct V2.0-anchored 704D residual expert. That native expert is experimental and not deployment-eligible. See [the local V2.0 integration record](docs/atom3d-v2-integration-local.md). Generic raw PDB/SDF-to-704D user inference is still not a portable, clean-room workflow.

## Architecture at a glance

- **Protein tower:** ESM-2 650M embeddings (1,280D) projected to 48D.
- **Compound tower:** 2,048-bit Morgan fingerprints (radius 2; chirality enabled) projected to 48D.
- **Feature fusion:** concatenation of protein, compound, elementwise product and absolute difference (192D).
- **Shared trunk:** 192 → 96 → 48, with task-specific scalar heads and frozen similarity blending.
- **Original public native structure expert:** separate 704D geometry residual on a legacy 423D V5 2D anchor; four source-trained folds. The local ATOM3D extension has a distinct 704D residual to its V2.0 ATOM3D head, but it is experimental and not deployment-eligible.

The released 2D inference implementation also contains a 128D auxiliary input named g, zero-filled by the portable Predictor; its training-time provenance is not fully documented. This is explicitly flagged rather than guessed. **GenSPARC is a benchmark/task head, not a GenSPARC input token.**

See [Architecture](docs/architecture.md) and [Methodology](docs/methodology.md).

## Headline results

### Double-cold comparisons

| Dataset | V2.0 RMSE | Matched reproduced PerceiverCPI RMSE |
|---|---:|---:|
| DAVIS (n=1,339) | **0.6216** | 0.6828 |
| KIBA (n=3,442) | **0.7864** | 0.9267 |
| METZ (n=1,416) | **0.8797** | 1.0187 |

These are comparisons against **locally reproduced** PerceiverCPI outputs, not the original paper's headline table values.

### GenSPARC comparison

V2.0 has **20/24 numerical metric wins** against rounded published RMSE/Pearson point estimates across 12 cohorts. These are not paired significance tests and do not establish universal superiority. Four metric misses are listed in [the benchmark ledger](results/benchmark-summary.md).

### Native geometry

- BioLiP3 pooled retrospective n=70: native RMSE 1.3184 vs calibrated 2D RMSE 1.5760.
- OpenBind single-target n=482: native RMSE 1.7986 vs calibrated 2D RMSE 1.7917 (slightly worse).
- The prespecified OpenBind REAL-vs-SHUFFLED MAE interval crosses zero.

**Do not claim universal 3D superiority or independent multi-target external validation.** The external validation gate is blocked for lack of an eligible new cohort.

## ATOM3D prior research and the V2.0 local extension

ATOM3D-LBA supported earlier GeoFuse structural R&D, including atom-pair descriptors, pocket and shuffled-geometry controls, scaffold-cold experiments and pose-sensitivity studies. Those V13?V36 experiments remain historical branches. The original public frozen V2.0 checkpoint did not include ATOM3D.

A separate local V2.0 extension has now been trained: it adds an ATOM3D pK head to the existing frozen ESM-2/Morgan shared representation. Its 704D native residual expert is retained for research only because it did not pass the geometry-specific deployment gate. The original DAVIS/KIBA/METZ/GenSPARC weights and headline results remain unchanged. The new checkpoint and ATOM3D-derived caches are local-only and are not part of the public release.

See [the historical ATOM3D research appendix](docs/atom3d-prior-research.md), [the local V2.0 integration record](docs/atom3d-v2-integration-local.md), the [historical aggregate summary](results/atom3d-historical-summary.json) and [evidence hashes](manifests/atom3d-evidence-manifest.json). No raw ATOM3D data, labels, structures, derived features, per-example predictions or ATOM3D-trained weights are redistributed.

To check the aggregate record without opening source data or rescoring the consumed test:

    python scripts/verify_atom3d_summary.py

## Reproducing the native analyses

The repository includes the frozen native expert implementation, label-free prepared feature/replay bundles and replay scripts. First verify exact prediction parity without accessing affinity labels:

    python reproduce.py --inference-only
    python reproduce_biolip3.py --inference-only

To recompute metrics and intervals from the original source affinity files:

    python -m pip install -r requirements-replay.txt
    python reproduce.py --archive data/OpenBind_EV-A71_2A.zip --out runs/openbind-reproduction.json
    python reproduce_biolip3.py --source data/experiments_affinity.tsv --out runs/biolip3-reproduction.json

The scripts verify source SHA-256 values and compare recomputed metrics/intervals with frozen references. See [Reproducibility](docs/reproducibility.md) and [Data provenance](docs/data-provenance.md).

## What is and is not independently reproducible

**Available now:** frozen checkpoint inference; raw 2D SMILES+FASTA path; prepared native expert inference; exact native prepared-prediction parity; OpenBind/BioLiP3 source-label replay scripts.

**Still not established:** original 2D training rerun from raw source data; full clean-room regeneration of ESM and native features; complete DAVIS/KIBA/METZ double-cold table regeneration from a fresh clone; GenSPARC raw-feature regeneration; independent clean-machine retraining. The auxiliary g training-time feature provenance also remains unresolved.

See [Reproducibility](docs/reproducibility.md) and [the reproduction gap register](docs/reproduction-gap-register.md).

## Repository boundaries and licensing

This repository contains the frozen V2.0 inference checkpoint and training-only anchor arrays, inference/replay code, label-free prepared replay inputs, aggregate results and methodology. It does not distribute raw affinity labels, raw benchmark structures, ESM-2 weights or the per-example PRED_NOY benchmark prediction archive.

OpenBind source data are CC0 1.0. BioLiP3 is freely available for non-commercial research; users must follow the source provider's terms. No license is assigned to this GeoFuse-CPI repository pending ownership and third-party compatibility review. Public visibility does not itself grant reuse rights. See [NOTICE](NOTICE.md) for source acknowledgements and data-use boundaries.
