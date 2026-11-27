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
6. The native prepared-feature replay is reproducible from the supplied label-free feature bundles. A separate predicted-structure CLI now generates a 704D descriptor from a docked pose, but that is not the same as reproducing the native expert's experimentally bound training-feature pipeline and is not applied to the native expert.

## Predicted-structure preparation limits

1. The optional raw sequence/SMILES route can obtain a sequence-matched AlphaFold DB model or, with explicit consent, submit a sequence to the public ESMFold API. These are predicted structures, not experimental structures.
2. RDKit conformers are starting ligand geometries. P2Rank predicts candidate pockets; its probability is not compound-specific binding probability.
3. AutoDock Vina produces approximate docking scores and poses. A Vina score is not pKd, Ki, Kd or experimental affinity.
4. The route can extract a 704D descriptor from a docked pose, but the descriptor is not fed into the current native residual expert. That expert was trained on experimentally bound complexes, failed its deployment gate, and has not been validated on docked poses.
5. No 3D affinity correction, docking accuracy claim, or 2D+3D superiority claim is made from this route. The 2D task head remains the affinity estimator.
6. Public ESMFold sequence submission is opt-in, rate-limited and unsuitable for bulk or confidential sequences. The current hosted-API path has a conservative 400-residue limit.
7. The CLI is a local research/inference tool; it is not yet wired into the production website frontend.
8. Docking uses a rigid receptor and the supplied SMILES protonation/tautomer state. It does not model protein flexibility, induced fit, solvent waters, metal/cofactor effects, or pH-dependent ligand-state enumeration.
9. P2Rank pocket ranking and Vina docking are heuristic candidate-generation stages; they require separate pose-quality and external affinity validation before being used for screening decisions.

## ATOM3D scope and test-set integrity

1. ATOM3D supported substantial earlier structural development. A separate local V2.0 extension now adds a fifth task-specific ATOM3D pK head; it does not replace the original public frozen checkpoint or alter the existing four benchmark heads.
2. The official 466-row validation and 490-row test were each scored once by earlier, non-V2.0 models. The one-time test model was a separate 3,343D LightGBM experiment.
3. The official test audit records that no reread or new scoring is permitted. Do not run V2.0 against that consumed test set or describe the earlier LightGBM score as a V2.0 result.
4. The native-versus-pocket and native-versus-shuffled confidence intervals in the one-time test crossed zero. The point-estimate ordering does not establish a geometry-specific benefit.
5. Earlier source-TRAIN experiments are retrospective and repeatedly used the same 3,507 training labels. Their confidence intervals are development diagnostics, not independent prospective evidence.
6. The historical evidence is documented in [atom3d-prior-research.md](atom3d-prior-research.md); the local V2.0 extension and its exact internal CV results are documented in [atom3d-v2-integration-local.md](atom3d-v2-integration-local.md).
7. The local ATOM3D head was evaluated only by retrospective, four-fold 3-mer-Jaccard protein-cluster CV within TRAIN. This is not independent validation or an aligned sequence-identity split. The native 704D residual failed its deployment gate because pocket-only was stronger; it must remain experimental.
8. The local ATOM3D-trained checkpoint and derived feature caches are not part of the public release. No raw ATOM3D data, labels, structures, derived feature arrays or per-example predictions are redistributed.

## External validation

The frozen gate is BLOCKED_NO_ELIGIBLE_COHORT. No new independent multi-target external cohort passed the predeclared eligibility and pose-quality requirements for scoring. Do not call existing retrospective BioLiP3/OpenBind results prospective external validation.

## Architecture and implementation limits

- V2.0 is a frozen hybrid release with a 2D protein/compound predictor and a separately trained native-structure expert. It is not one jointly trained 2D+3D network.
- The 128-dimensional auxiliary input g in the 2D inference class is zero-filled by the released Predictor, and its training-time semantics are not fully documented.
- The original complete V2.0 2D training script and full training pipeline were not present in the preserved release checkout.
- The public repository includes the frozen checkpoint and training-only anchors, but not raw affinity datasets, protected labels, or ESM-2 weights.
- The optional predicted-structure CLI can generate a docked-pose 704D descriptor, but it does not provide a validated 3D affinity prediction and is not yet wired into the production website frontend.
- Full raw-input portability for all GenSPARC/native feature paths is not established.

## Reproducibility interpretation

The release now supports raw 2D inference, prepared-feature native inference, exact native prediction parity, and source-label replay for OpenBind and BioLiP3. These are meaningful reproducibility improvements, but they are not independent retraining. The DAVIS/KIBA/METZ double-cold tables and GenSPARC results cannot yet be regenerated end-to-end from a fresh clone.

## Version isolation

Only V2.0 results belong in this repository's headline ledger. V363, V499, V48, v6, later fusion experiments and website production results are separate versions and must not be score-stitched into V2.0.
