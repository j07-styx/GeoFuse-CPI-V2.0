# V2.0 benchmark summary

All values below are from the single frozen V2.0 checkpoint identified in [the artifact identity record](../manifests/frozen-artifact-identity.json). See the machine-readable records in this directory.

## DAVIS, KIBA and METZ double-cold

Lower RMSE is better. The comparator is a matched local PerceiverCPI reproduction, not the original published PerceiverCPI table value.

| Dataset | n | V2.0 RMSE | Reproduced PerceiverCPI RMSE | Difference (V2.0 − comparator) |
|---|---:|---:|---:|---:|
| DAVIS | 1,339 | **0.6216** | 0.6828 | −0.0611 |
| KIBA | 3,442 | **0.7864** | 0.9267 | −0.1402 |
| METZ | 1,416 | **0.8797** | 1.0187 | −0.1390 |

These are favorable matched-reproduction point estimates. They are **not** verified wins over the original PerceiverCPI published headline values.

## GenSPARC published-point-estimate comparison

V2.0 has 20/24 numerical metric wins across 12 cohorts × 2 metrics (RMSE and Pearson), against the paper's rounded published values. Four metrics are losses:

| Cohort | Metric | V2.0 | Published |
|---|---|---:|---:|
| Original unseen compound | Pearson | 0.7835 | 0.787 |
| Original unseen both | RMSE | 1.5527 | 1.552 |
| Sequence-hard unseen protein | RMSE | 1.5465 | 1.444 |
| Sequence-hard unseen protein | Pearson | 0.5824 | 0.642 |

These are not paired significance tests. The 20/24 count is not 20/24 independent datasets and is not universal superiority over all methods in the GenSPARC paper.

## Native geometry ablations

| Cohort | n | Calibrated 2D RMSE | Native RMSE | Shuffled RMSE | Pocket-only RMSE |
|---|---:|---:|---:|---:|---:|
| BioLiP3 first | 37 | 1.5209 | **1.2232** | 1.2631 | 1.3158 |
| BioLiP3 second | 33 eligible | 1.6356 | **1.4176** | 1.4289 | 1.4523 |
| BioLiP3 pooled post hoc | 70 | 1.5760 | **1.3184** | 1.3438 | 1.3818 |
| OpenBind single-target | 482 | **1.7917** | 1.7986 | 1.8222 | 1.9162 |

BioLiP3 shows favorable retrospective native-vs-2D point estimates. OpenBind native geometry is slightly worse than calibrated 2D on RMSE. The native-vs-shuffled primary MAE criterion on OpenBind did not pass its prespecified interval test. Do not claim universal 3D superiority.

A separate raw v5 2D arm on OpenBind had RMSE 0.8803; it is retained as a distinct raw arm in the audit and should not be silently substituted for the calibrated 2D comparator in the table above.

## External validation

The frozen external gate is `BLOCKED_NO_ELIGIBLE_COHORT`. No independent new multi-target external cohort is authorized for scoring under the preregistered gate. Current BioLiP3/OpenBind analyses are retrospective and do not satisfy that gate.

## Interpretation

V2.0 has promising multi-dataset 2D benchmark evidence and useful native-structure ablations. It has not established independent external superiority, universal geometry-specific gains, or an end-to-end jointly trained 2D+3D architecture.
