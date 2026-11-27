# Release scope

## Included in this repository

- V2.0 architecture, methodology, aggregate benchmark ledger and limitations.
- Frozen checkpoint GEOFUSE_CPI_V2.pt and label-free training-only anchors TRAIN_ANCHORS_NO_Y.npz, distributed inside the checksum-verified inference asset archive.
- Raw SMILES + FASTA 2D inference implementation using ESM-2 650M.
- Prepared-feature native expert inference implementation.
- ESM-2 download helper with a pinned SHA-256.
- OpenBind and BioLiP3 label-free prepared-feature replay bundles and exact inference-parity scripts.
- Source-label replay scripts, aggregate comparator references and environment pins.
- Public metadata checker and frozen-inference verification script.
- Historical ATOM3D appendix and aggregate-only evidence manifest.
- ATOM3D V2.0 training/inference code and a local integration record. The trained ATOM3D checkpoint and all ATOM3D-derived feature caches remain excluded from the public release.

## Excluded

- Raw affinity labels and protected assay outcomes.
- Raw OpenBind/BioLiP3 PDB/SDF structures and original source datasets.
- Raw ATOM3D source archives, affinity labels, structures, derived feature arrays, per-example predictions and ATOM3D-trained checkpoint weights. The local research checkpoint is ignored by Git and not redistributed.
- ESM-2 650M weights (downloaded from the official fair-esm host).
- The PRED_NOY.npz per-example benchmark prediction archive.
- The original complete 2D training script and full training data/split pipeline, which were not present in the preserved V2.0 release checkout.
- Private local virtual environments, credentials, machine-specific paths and unrelated V363/V499/V48/v6 experimental artifacts.

## Intended use

This is a frozen research inference and evidence release. It is not the production website backend or a hosted service. The raw 2D predictor accepts new valid SMILES and protein FASTA inputs. The native expert accepts prepared 423D/704D features; it does not create poses or infer native geometry from arbitrary raw structures.

## Licensing

No license is assigned to this repository at this time. OpenBind source data are CC0 1.0. BioLiP3 source data are freely available for non-commercial research use; users must comply with the provider's terms. ESM-2 weights are not redistributed. Public visibility does not itself grant a license to reuse this repository's code or model weights.
