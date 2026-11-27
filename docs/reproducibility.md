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

### 6. ATOM3D V2.0 local extension boundary

A local V2.0 extension adds an ATOM3D pK head on the frozen ESM-2/Morgan representation and a 704D native residual expert. It was trained and evaluated only on the official 3,507-row TRAIN split using nested four-fold heuristic protein-cluster CV. The official validation/test were not opened; the test was already consumed by an earlier non-V2.0 model and must not be rescored. The native expert failed the deployment gate; pocket-only was stronger. This is retrospective development evidence, not independent validation or a headline V2.0 benchmark. See [the local integration record](atom3d-v2-integration-local.md).

### 7. Full benchmark retraining and clean-room reproduction

**Not established.** The repository does not contain the original 2D training script, all original raw source feature files, complete cross-source split manifests, or all raw benchmark inputs needed to recreate every DAVIS/KIBA/METZ/GenSPARC table from scratch. The original local release verifier and smoke test passed on the research machine, but this is not equivalent to an independent clean-machine retraining run.

### 8. Predicted-structure preparation and docking smoke tests

The optional structure-preparation CLI is implemented in `inference/structure_pipeline.py`. The local Windows environment used P2Rank 2.5.1, AutoDock Vina 1.2.7, Meeko 0.8.0, Gemmi 0.7.5, RDKit 2025.09.2, Biopython 1.88 and SciPy 1.17.1. P2Rank and Vina binaries are downloaded by `scripts/install_structure_tools.ps1`, SHA-256 checked, and kept outside Git.

The 2D prediction path was run with aspirin and `examples/AAK1.fasta` using both checkpoints and the same cached, checksum-pinned ESM-2 650M weights:

- Public frozen V2.0 DAVIS output: 5.476265008935864 pKd.
- Local ATOM3D-extended checkpoint DAVIS output: 5.476265008935864 pKd.
- Absolute difference: exactly 0.0. This confirms the existing DAVIS head output was preserved for this smoke input; it is not an accuracy or affinity-measurement claim.

The structure route was smoke-tested on the public 142-residue human hemoglobin-alpha sequence with aspirin as a **software integration test only**:

| Protein structure source | Sequence identity / coverage | Mean pLDDT | Top P2Rank probability | Vina score | Geometry |
|---|---:|---:|---:|---:|---|
| AlphaFold DB API (P69905) | 1.00 / 1.00 | 98.06 | 0.964 | -5.947 kcal/mol | finite 704D |
| ESMFold API (explicit opt-in) | 1.00 / 1.00 | 94.69 | 0.964 | -5.998 kcal/mol | finite 704D |

A separate user-supplied AlphaFold PDB test also completed and generated a finite 704D descriptor. The tests demonstrate that structure acquisition, sequence matching, P2Rank, Meeko preparation, Vina docking, pose reconstruction and descriptor extraction execute end-to-end. They do **not** demonstrate that aspirin binds hemoglobin, that the predicted poses are accurate, or that Vina scores correlate with measured affinity.

The current ATOM3D 704D native residual expert is deliberately **not applied** to these docked poses. It was trained on experimentally bound complexes, failed its deployment gate against the pocket-only control, and has not been validated on predicted-pose inputs. The pipeline therefore returns a 2D task-head estimate plus a computationally predicted complex and 704D research descriptor—not a 3D-corrected affinity.

The ESMFold API receives the protein sequence only when the user passes `--allow-public-esmfold`. The public service is rate-limited and should not be used for bulk folding or confidential sequences. The current hosted route has a conservative 400-residue cap. A local PDB/mmCIF or a sequence-matched AlphaFold DB UniProt accession avoids sending the raw sequence to ESMFold.

See [structure-aware-inference.md](structure-aware-inference.md) for installation, commands and output descriptions.

## Important architecture/reproducibility caveat: auxiliary g

The portable 2D model class has a 128-dimensional auxiliary input g; the published raw-input wrapper supplies zeros. The checkpoint contains nonzero weights for the g projection, but the original training-time provenance and values for this input are not recoverable from the public inference bundle. Do not invent a biological interpretation for it, and do not claim the full training pipeline is reproduced until the original source is recovered and audited.

## Comparator boundaries

- DAVIS/KIBA/METZ results compare against matched locally reproduced PerceiverCPI outputs, not the original PerceiverCPI paper's headline values.
- GenSPARC comparisons are against rounded published point estimates, not paired per-example predictions.
- Native geometry is a separate residual expert on a legacy V5 anchor, not an end-to-end jointly trained continuation of the V2.0 2D trunk.
- BioLiP3 and OpenBind results are retrospective. The independent external validation gate remains blocked.
