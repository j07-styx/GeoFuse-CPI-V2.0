# Reproducibility and verification

## Scope

This repository contains the frozen V2.0 checkpoint and training-only similarity anchors, raw-input 2D inference code, prepared-feature native inference code, label-free replay bundles, source-label replay scripts and replay reports. It is still **not** a complete clean-room retraining release.

## Frozen artifact identity

- Checkpoint SHA-256: b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4
- Training-anchor SHA-256: 42ddd4d3af1f00f15bf96985adaaa15b821fb67fe6a5f2fddf3fa796a4773769
- Inference asset archive SHA-256: 71148521ae43187fcce4183cd7d388dc79daa04ebd6ace98bfacc08685a37017
- ESM-2 650M weight SHA-256: ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c

The public inference archive contains only GEOFUSE_CPI_V2.pt and TRAIN_ANCHORS_NO_Y.npz. It deliberately excludes PRED_NOY.npz, raw labels and source assay tables.

## Reproduction levels

### 1. Public-record consistency

From the repository root:

    python scripts/verify_public_records.py

This checks aggregate-record and artifact-identity metadata. It does not reproduce predictions.

### 2. Frozen raw-input 2D inference

Install the pinned reference dependencies and prepare assets:

    python -m pip install -r inference/requirements-inference.txt
    python scripts/prepare_inference_assets.py
    python scripts/download_esm2_weights.py

Run the bundled raw-input smoke test:

    python scripts/verify_frozen_inference.py --raw

The smoke example uses aspirin and the AAK1 FASTA sequence. On the original Windows/Python 3.11.9 CPU environment, the frozen raw-input predictor returned DAVIS pKd 5.4762650089. The stored value is a software parity fixture, **not** a measured affinity or accuracy claim.

The raw predictor uses fair-esm 2.0.0 and the public ESM-2 650M checkpoint. The ESM weights are downloaded from the official fair-esm model host and SHA-256 checked before loading.

### 3. Prepared native expert inference

The same frozen checkpoint contains four native residual folds. Run:

    python -m pip install -r inference/requirements-inference.txt
    python -m pip install -r inference/requirements-native.txt
    python scripts/prepare_inference_assets.py
    python inference/predict_native_prepared.py --features path/to/prepared_noY.npz --out predictions.csv

Input arrays must be label-free and have exact shapes x2d=(N,423) and real/shuffled/pocket=(N,704). The script rejects fields named y, label, affinity or target.

This is prepared-feature inference only. It does not generate native ligand poses or derive the 423D/704D descriptors from raw PDB/SDF files.

### 4. Exact native prepared-prediction parity

The OpenBind and BioLiP3 replay bundles contain label-free prepared features and frozen reference predictions. They do not contain numerical affinity labels or raw PDB/SDF coordinates.

    python reproduce.py --inference-only
    python reproduce_biolip3.py --inference-only

Verified on the original Windows/Python 3.11.9 CPU environment and re-run from the release branch:

- OpenBind: 660 prepared inputs; maximum absolute difference was exactly 0.0 for 2D, native, shuffled, pocket-only and raw V5 arms.
- BioLiP3: 73 prepared inputs; maximum absolute difference was exactly 0.0 for the same arms.

The native fold weights and calibrators in the OpenBind replay bundle were also compared directly with the frozen checkpoint; the tested fold tensors and normalization arrays matched exactly.

### 5. Source-label metric replay

**Completed on the original Windows/Python 3.11.9 CPU environment and repeated from the release branch.** The scripts downloaded the original source files, verified their pinned hashes, recomputed the metrics and confidence intervals, and matched the frozen aggregate references.

OpenBind:
- Source: Zenodo record 22142262, 124,498,453 bytes, SHA-256 e0d2e38eddf75ebf4247aff0fe28ae15dfa25ecb30deeb4cf950ae184b304a25.
- Scored cohort: 482 compounds across 264 scaffold clusters.
- Verification: 60 metric entries and 12 bootstrap intervals passed.
- Calibrated 2D RMSE: 1.791742.
- Native REAL geometry RMSE: 1.798589.
- SHUFFLED geometry RMSE: 1.822189.
- POCKET-only RMSE: 1.916216.
- Native REAL is slightly worse than calibrated 2D on RMSE. The prespecified REAL-minus-SHUFFLED MAE interval crosses zero.

BioLiP3:
- Source: official BioLiP experiments_affinity.tsv, 36,087,375 bytes, SHA-256 c8467ed78b4e78889b417d752e8da4973bfd3a95fb5dc6f1244cafcc8249d0ab.
- Cohorts: 37 and 33 eligible measured cases; 70 pooled post hoc.
- Verification: 105 metric entries and 8 stratified confidence intervals passed.
- Pooled native REAL RMSE: 1.318425.
- Pooled SHUFFLED geometry RMSE: 1.343804.
- Pooled calibrated 2D RMSE: 1.576021.
- REAL-minus-SHUFFLED MAE 95% interval: [-0.0564, 0.0248], which crosses zero.

Machine-readable run reports are in results/reproduction/ and the source/bundle hashes and environment are recorded in reproduction/replay-run-manifest.json.

These are retrospective source-label replays from frozen prepared features and fixed model artifacts. They are not independent retraining, raw feature regeneration, prospective validation or proof of universal 3D superiority.

### 6. ATOM3D historical result boundary

ATOM3D is not a verified V2.0 evaluation cohort. The official validation and test sets were each scored once by earlier non-V2.0 models; the test was scored with a separate 3,343D LightGBM model. Its final integrity audit prohibits rereading or new scoring. The archived point estimates and confidence intervals are documented in [the ATOM3D appendix](atom3d-prior-research.md), but they are not V2.0 scores and cannot be used as independent validation of the V2.0 checkpoint.

### 7. Full benchmark retraining and clean-room reproduction

**Not established.** The repository does not contain the original 2D training script, all original raw source feature files, complete cross-source split manifests, or all raw benchmark inputs needed to recreate every DAVIS/KIBA/METZ/GenSPARC table from scratch. The original local release verifier and smoke test passed on the research machine, but this is not equivalent to an independent clean-machine retraining run.

## Important architecture/reproducibility caveat: auxiliary g

The portable 2D model class has a 128-dimensional auxiliary input g; the published raw-input wrapper supplies zeros. The checkpoint contains nonzero weights for the g projection, but the original training-time provenance and values for this input are not recoverable from the public inference bundle. Do not invent a biological interpretation for it, and do not claim the full training pipeline is reproduced until the original source is recovered and audited.

## Comparator boundaries

- DAVIS/KIBA/METZ results compare against matched locally reproduced PerceiverCPI outputs, not the original PerceiverCPI paper's headline values.
- GenSPARC comparisons are against rounded published point estimates, not paired per-example predictions.
- Native geometry is a separate residual expert on a legacy V5 anchor, not an end-to-end jointly trained continuation of the V2.0 2D trunk.
- BioLiP3 and OpenBind results are retrospective. The independent external validation gate remains blocked.
