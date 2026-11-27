# Native replay bundles

These archives are prepared-feature replay inputs, not raw structure datasets.

## OpenBind

- Archive: openbind-replay.zip
- SHA-256: 3dc945c3e68d9459319460f32115ede09688fb1b692b97c4f12baba3c09d9984
- Contains 660 label-free feature rows, four frozen native-expert folds and frozen reference predictions.
- Source dataset: https://doi.org/10.5281/zenodo.22142262 (CC0 1.0).

## BioLiP3

- Archive: biolip3-replay.zip
- SHA-256: 85aecd9a195e5fa93046bfe20be49aa324cd673a9d14ee25a06f460d963ee607
- Contains two prepared cohorts (37 + 36 rows), label-free feature arrays and frozen reference predictions. The scoring script selects 37 + 33 measured cases from the source affinity table.
- Source data: https://zhanggroup.org/BioLiP/lh_3affinity_data/00_Experiments/experiments_affinity.tsv.
- BioLiP3 is freely available for non-commercial research use. Follow the source provider's terms; this repository does not grant additional rights.

## Exact prediction parity, without labels

From the repository root, after installing requirements-replay.txt:

    python reproduce.py --inference-only
    python reproduce_biolip3.py --inference-only

Expected results on the pinned reference environment:

- OpenBind: 660 rows; maximum absolute difference 0.0 for each frozen prediction arm.
- BioLiP3: 73 rows; maximum absolute difference 0.0 for each frozen prediction arm.

## Recompute source-label metrics

The scripts download or accept the source files, verify their pinned SHA-256 values, compute metrics and confidence intervals, and compare the outputs to the frozen aggregate references.

    python reproduce.py --archive data/OpenBind_EV-A71_2A.zip --out runs/openbind-reproduction.json
    python reproduce_biolip3.py --source data/experiments_affinity.tsv --out runs/biolip3-reproduction.json

The replay is retrospective and structure-conditioned. It is not prospective pose-free screening and does not establish universal 3D superiority.
