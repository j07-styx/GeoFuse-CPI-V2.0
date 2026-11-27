# Reproduction gap register

This register distinguishes completed release work from unresolved research reproducibility requirements.

| Item | Status | Evidence / remaining action |
|---|---|---|
| Frozen V2.0 checkpoint publicly accessible | Complete | Checksum-pinned inference asset archive contains the checkpoint. |
| Training-only similarity anchors accessible | Complete | Checksum-pinned, label-free TRAIN_ANCHORS_NO_Y.npz included. |
| Raw SMILES + FASTA 2D inference | Complete on reference machine | ESM-2 650M raw-input smoke output reproduced to tolerance 0.001. |
| Prepared native expert inference | Complete on reference machine | Frozen fourfold native expert exposed through a prepared-feature CLI. |
| Native replay prediction parity | Complete on reference machine | 660 OpenBind and 73 BioLiP3 rows; maximum absolute prediction difference 0.0. |
| OpenBind source-label metric replay | Complete on reference machine and release branch | 482 compounds; 60 metric checks and 12 bootstrap intervals passed against frozen aggregates. |
| BioLiP3 source-label metric replay | Complete on reference machine and release branch | 37 + 33 cohorts; 105 metric checks and 8 stratified intervals passed against frozen aggregates. |
| DAVIS/KIBA/METZ complete double-cold replay | Incomplete | Need exact raw source files, feature generation code, row-level split manifests and evaluation inputs for the frozen V2.0 run. |
| GenSPARC full raw-feature replay | Incomplete | Published point-estimate comparisons are documented; original per-example inputs and raw feature-generation pipeline are not included. |
| Original V2.0 2D training script | Not recovered | The preserved release has checkpoint/inference code but not the original complete training program. Do not substitute a newly written trainer and call it exact reproduction. |
| Auxiliary g training-time provenance | Unresolved | Portable inference zero-fills g; original training-time feature source/values are not in the release manifest. |
| Native raw PDB/SDF-to-704D feature pipeline | Incomplete | Native model inference is available for prepared features; end-to-end portable feature generation still needs packaging and validation. |
| Independent clean-machine retraining | Not established | Reference-machine smoke tests and source-label replay are not a clean-room retraining study. |

## Release claim

The accurate claim is: **GeoFuse-CPI V2.0 has a checksum-verified frozen inference release with runnable raw 2D inference and prepared-feature native inference. Its BioLiP3 and OpenBind native results have exact prepared-prediction parity and source-label metric replay scripts that pass against frozen references. It is not yet a fully independently retrainable release, and not every headline benchmark table can be regenerated from a fresh clone alone.**
