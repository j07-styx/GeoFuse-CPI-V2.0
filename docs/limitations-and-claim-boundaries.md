# Limitations and claim boundaries

These constraints are part of the research record and must accompany any presentation or manuscript using V2.0 results.

## Benchmark comparison limits

1. DAVIS/KIBA/METZ comparisons are against matched locally reproduced PerceiverCPI outputs, **not** original PerceiverCPI published headline table values.
2. GenSPARC comparisons are against rounded published point estimates; per-interaction predictions were unavailable, so paired significance testing was not possible.
3. The 20/24 GenSPARC metric count is not 20/24 independent datasets and does not establish universal superiority over all methods in that paper.
4. Cross-source training differs from original published benchmark protocols. Exact entity quarantine is not full scaffold or protein-homology isolation. Historical cohorts were inspected during research.

## Native geometry limits

1. The native expert is separately trained and uses a legacy 423-dimensional V5 anchor. It is not end-to-end jointly trained with the V2.0 2D trunk.
2. BioLiP3 is retrospective, mixes Ki/Kd, and its pooled 70 is post hoc. The native-vs-shuffled geometry-specific effect is not established as universal.
3. OpenBind is one protease family. Native geometry RMSE (1.7986) is slightly worse than calibrated 2D RMSE (1.7917).
4. On OpenBind, the prespecified REAL-vs-SHUFFLED MAE interval crosses zero. A favorable ranking comparison does not rescue the failed primary absolute-error criterion.
5. Experimentally bound native poses make this a structure-conditioned rescoring experiment, not pose-free prospective screening.
6. The native prepared-feature replay is reproducible from the supplied label-free feature bundles. This does not establish a portable raw PDB/SDF-to-704D feature-generation workflow.

## External validation

The frozen gate is BLOCKED_NO_ELIGIBLE_COHORT. No new independent multi-target external cohort passed the predeclared eligibility and pose-quality requirements for scoring. Do not call existing retrospective BioLiP3/OpenBind results prospective external validation.

## Architecture and implementation limits

- V2.0 is a frozen hybrid release with a 2D protein/compound predictor and a separately trained native-structure expert. It is not one jointly trained 2D+3D network.
- The 128-dimensional auxiliary input g in the 2D inference class is zero-filled by the released Predictor, and its training-time semantics are not fully documented.
- The original complete V2.0 2D training script and full training pipeline were not present in the preserved release checkout.
- The public repository includes the frozen checkpoint and training-only anchors, but not raw affinity datasets, protected labels, or ESM-2 weights.
- Full raw-input portability for all GenSPARC/native feature paths is not established.

## Reproducibility interpretation

The release now supports raw 2D inference, prepared-feature native inference, exact native prediction parity, and source-label replay for OpenBind and BioLiP3. These are meaningful reproducibility improvements, but they are not independent retraining. The DAVIS/KIBA/METZ double-cold tables and GenSPARC results cannot yet be regenerated end-to-end from a fresh clone.

## Version isolation

Only V2.0 results belong in this repository's headline ledger. V363, V499, V48, v6, later fusion experiments and website production results are separate versions and must not be score-stitched into V2.0.
