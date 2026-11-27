# ATOM3D integration in GeoFuse-CPI V2.0

**Status:** Local V2.0 research checkpoint created. This is not a new model generation, not V4, and not the public frozen V2.0 release artifact.

## What was integrated

The original V2.0 checkpoint remains unchanged. A local extension was created at:

`models/v2.0/GEOFUSE_CPI_V2_ATOM3D_LOCAL.pt`

The extension adds an ATOM3D-specific affinity head to the existing V2.0 shared representation:

- Frozen ESM-2 650M protein tower.
- Existing Morgan 2048-bit compound tower.
- Existing concatenation, elementwise-product and absolute-difference fusion.
- New `atom3d` pK task head, fitted as a group-balanced Ridge head over the frozen 48D fused representation.
- Optional 704D native ligand-pocket residual expert, stored in the same checkpoint but marked experimental.

The original DAVIS, KIBA, METZ and GenSPARC tensors were verified bit-for-bit unchanged. The original `GEOFUSE_CPI_V2.pt` remains in place and correctly refuses an `atom3d` prediction because it has no ATOM3D head.

This is a fifth task-specific head within V2.0, not a change to the shared trunk. That choice preserves the published benchmark heads and their existing scores. Fine-tuning the shared trunk on ATOM3D would require retraining and re-evaluating all existing task heads.

## ATOM3D data used

The ATOM3D LBA 30%-sequence-identity split contains 3,507 training, 466 validation and 490 test complexes. The task target is pK = -log10(K), with Ki/Kd affinity records. See the [official LBA description](https://www.atom3d.ai/lba.html) and [Zenodo dataset record](https://doi.org/10.5281/zenodo.4914718).

For this V2.0 integration:
- All 3,507 official TRAIN records were prepared.
- 2,105 unique observed protein sequences were encoded with the same ESM-2 650M layer-33 embedding convention used by V2.0.
- All 3,507 compounds received Morgan fingerprints: 2,961 from standard sanitized RDKit parsing and 546 from an explicitly flagged unsanitized-valence fallback. No source compounds were silently discarded.
- Each complex received a 704D native ligand-pocket descriptor plus shuffled-geometry and pocket-only controls.
- A four-fold split was made from 1,004 connected components using observed-chain 3-mer Jaccard >= 0.28. This is a heuristic sequence-cluster grouping, not a true aligned sequence-identity split.

**Official ATOM3D validation and test were not opened or rescored during this integration.** The official test had already been scored in an earlier, separate non-V2.0 experiment; it must not be reused as a fresh blind test.

The ATOM3D TRAIN labels were used in earlier V13–V36 exploratory work. The current cross-validation is therefore retrospective development evidence, not independent prospective validation.

## Selected ATOM3D 2D head

The group-balanced head was selected because it improved the heuristic protein-group macro metric compared with an unweighted head, with a small pooled-RMSE trade-off.

| Head, nested four-fold TRAIN CV | Pooled RMSE | Pearson | Protein-group macro RMSE |
|---|---:|---:|---:|
| Unweighted V2.0 ATOM3D head | 1.5578 | 0.6290 | 1.1850 |
| **Group-balanced V2.0 ATOM3D head** | **1.5666** | **0.6324** | **1.1713** |

The selected head is trained on the full 3,507-complex TRAIN source after the nested evaluation. These metrics are not comparable to the official ATOM3D validation/test leaderboard because they use an internal heuristic group-CV protocol.

## Native-geometry evaluation

The 704D native residual expert was trained with group-balanced loss and nested out-of-fold 2D residual targets. The same outer folds were used for all arms.

| Arm | Pooled RMSE | Pearson | Protein-group macro RMSE |
|---|---:|---:|---:|
| V2.0 2D group-balanced | 1.5666 | 0.6324 | **1.1713** |
| 2D + native 704D geometry | 1.4974 | 0.6776 | 1.1975 |
| 2D + shuffled 704D geometry | 1.5359 | 0.6619 | 1.2157 |
| 2D + pocket-only control | **1.4775** | **0.6815** | **1.1491** |

Paired group-bootstrap 95% intervals for native-minus-comparator group-macro RMSE:

- Native vs 2D: -0.0074 to 0.0602 (crosses zero).
- Native vs shuffled: -0.0510 to 0.0143 (crosses zero).
- Native vs pocket-only: 0.0207 to 0.0763 (native is worse).

The native expert therefore **fails the prespecified deployment gate**. The pocket-only control performs better on both pooled and group-macro RMSE. The native expert is retained in the checkpoint for research diagnostics only and is marked `deployment_eligible=false`.

Additional group-balanced controls also failed to establish a native-specific benefit:
- 1,210D residue-aware native geometry: pooled RMSE 1.5214, group-macro RMSE 1.2218.
- 55D Vina-like physical descriptors: pooled RMSE 1.5728, group-macro RMSE 1.1811.
- Combined 1,265D residue-aware plus physical descriptors: pooled RMSE 1.5068, group-macro RMSE 1.2168.

The 55D descriptor is Vina-like, not a calibrated AutoDock Vina energy. Pocket-only controls use the experimental contact graph and are diagnostic controls, not deployable pose-free features.

## Inference

The new checkpoint supports the ATOM3D 2D task head:

```python
import sys
sys.path.insert(0, "inference")
from geofuse_portable_kinase import Predictor

predictor = Predictor(
    "models/v2.0/GEOFUSE_CPI_V2_ATOM3D_LOCAL.pt",
    "models/v2.0/TRAIN_ANCHORS_NO_Y.npz",
)
pK = predictor.predict("atom3d", smiles, esm2_embeddings)
```

The native-geometry helper is `inference/predict_atom3d_native_prepared.py`. It requires a precomputed 704D descriptor from an experimentally bound ligand-pocket complex. It does not generate a pose or perform docking. Because the expert failed its validation gate, it refuses normal inference unless `--allow-experimental` is explicitly supplied.

## Separate predicted-structure preparation route

The V2.0 repository now also includes an optional raw sequence/SMILES preparation workflow in `inference/structure_pipeline.py`. It can retrieve a sequence-matched AlphaFold DB model, use an explicitly supplied PDB/mmCIF, or opt in to the public ESMFold API; then it generates RDKit ligand conformers, predicts pockets with P2Rank, docks with AutoDock Vina, and extracts a 704D descriptor from the selected docked pose.

This is a separate inference-preparation workflow, not a new ATOM3D training cohort or a change to the V2.0 shared architecture. The generated docked-pose descriptor is **not** fed into the ATOM3D native residual expert and no 3D affinity correction is returned. The native expert was trained on experimentally bound complexes and failed its deployment gate; it is not validated on predicted/docked poses. Full details and software smoke-test boundaries are in [structure-aware-inference.md](structure-aware-inference.md).

## Claim boundaries

- ATOM3D is now represented by an additional task-specific head in a local GeoFuse-CPI V2.0 checkpoint.
- The original V2.0 DAVIS/KIBA/METZ/GenSPARC heads and their weights are unchanged.
- No official ATOM3D validation or test result is claimed for this checkpoint.
- No PerceiverCPI head-to-head win is claimed from this internal CV.
- No native-3D superiority is claimed; pocket-only was stronger.
- The internal folds are heuristic 3-mer sequence clusters and are not a replacement for independent external validation.
- The historical official test remains consumed and must not be reopened.

## Licensing and distribution

The ATOM3D LBA dataset is marked CC BY-NC-ND by its official project page and Zenodo record. The ATOM3D-trained checkpoint and feature caches are kept local and are not committed or uploaded. Do not redistribute the trained weights or use them commercially without permission and an appropriate rights review.

## Local artifact identity

- Checkpoint: `models/v2.0/GEOFUSE_CPI_V2_ATOM3D_LOCAL.pt`
- SHA-256: `de5e78d40388f001c6e998af487f3146c5742954b20e1d5eacbe1926224527ea`
- Size: 7,287,411 bytes
- The original public baseline remains `models/v2.0/GEOFUSE_CPI_V2.pt`.
