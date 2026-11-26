# Methodology and evaluation protocol

## Frozen candidate

All V2.0 headline values refer to one frozen hybrid research checkpoint:

- Checkpoint SHA-256: `b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4`
- Prediction archive SHA-256: `b8f935e2115ead2b83e247477118adedd38559acf368b80eabafc160240f0555`

The results ledger explicitly excludes later experimental variants. Do not combine V363, V499, v6, V48 or any other version's scores with V2.0.

## 2D training/evaluation cohorts

| Dataset | Rows before filtering | Rows after filtering | Validation rows | Held-out evaluation |
|---|---:|---:|---:|---:|
| DAVIS | 12,780 | 9,144 | 550 | Double-cold: 1,339 |
| KIBA | 53,045 | 27,299 | 3,166 | Double-cold: 3,442 |
| METZ | 15,320 | 7,288 | 658 | Double-cold: 1,416 |
| GenSPARC | 2,104 | 2,017 | 233 | 12 reported cohorts |

The 12 GenSPARC cohorts comprise original, sequence-hard and structure-hard partitions, each split into seen-both, unseen-compound, unseen-protein and unseen-both. Their sizes are recorded in `results/PAPER_RESULTS.json`.

The cross-source training design differs from the original published benchmark protocols. The frozen record describes exact entity quarantine; it does **not** establish full scaffold isolation or protein-homology isolation across every source. Historical benchmark cohorts were inspected during model development. Treat the current results as retrospective evidence, not a new untouched prospective test.

## 2D representations and predictor

- Protein sequence representation: frozen ESM-2 650M, 1,280-dimensional embeddings.
- Compound representation: 2,048-bit Morgan fingerprints, radius 2, chirality enabled.
- Learned protein and compound projections, interaction feature concatenation, shared trunk and assay-specific heads: see [architecture.md](architecture.md).
- The inference implementation blends the neural output with a task-specific similarity estimate using frozen training anchors and kernel weights.

## Native-structure evaluation

The native expert uses a legacy 423-dimensional 2D anchor and 704-dimensional geometry descriptors. It is trained separately from the V2.0 2D trunk. Four source-trained folds are used. The reported controls are calibrated 2D anchor, native geometry, shuffled geometry and pocket-only geometry.

BioLiP3 reports two retrospective cohorts with 37 and 33 eligible measured cases (70 pooled post hoc). These cohorts are not independent assay databases; the pooled score is computed from pooled MSE, not by averaging cohort RMSEs.

OpenBind is a single-target EV-A71/CVA16 2A protease source. 660 source-only compounds passed structural/chemistry selection and were predicted before numerical pKd was opened; 482 had measured pKd values for scoring. It is not a multi-target external benchmark.

## Comparator rules

- DAVIS/KIBA/METZ comparisons are against matched **locally reproduced** PerceiverCPI outputs, not the original PerceiverCPI paper's headline table values.
- GenSPARC comparisons are against rounded published point estimates. Per-interaction published predictions were unavailable, so these are point-estimate comparisons, not paired significance tests.
- Native 3D results are retrospective structure-conditioned rescoring results. They must not be presented as prospective pose-free screening performance.

## Model-selection rule

The V2.0 paper-results record states that one existing frozen candidate was selected before later fusion experiments. No test-specific checkpoint selection or score stitching is permitted. Any replacement checkpoint must be re-evaluated across all headline cohorts before its results can replace this ledger.
