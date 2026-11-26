# Dataset and source provenance

## Affinity benchmark sources

The V2.0 result ledger includes DAVIS, KIBA, METZ and GenSPARC. Exact row counts, filtering counts, validation sizes and held-out cohort sizes are recorded in [PAPER_RESULTS.json](../results/PAPER_RESULTS.json).

- DAVIS, KIBA and METZ are affinity benchmark sources with different endpoint conventions. KIBA's composite score is not interchangeable with pKd.
- GenSPARC source paper: [Communications Chemistry (2025)](https://doi.org/10.1038/s42004-025-01844-0). The V2.0 comparison uses rounded published RMSE/Pearson point estimates.
- Original PerceiverCPI source implementation: [dmis-lab/PerceiverCPI](https://github.com/dmis-lab/PerceiverCPI), revision `755f076ef51099a90332fd1ec8437f94e3750c0b`.

## Protein and compound features

- Protein sequence embeddings: ESM-2 650M, 1,280 dimensions.
- Compound fingerprints: RDKit Morgan fingerprints, radius 2, 2,048 bits, chirality enabled.
- The frozen predictor also contains an auxiliary 128-dimensional input `g`; the released portable inference wrapper supplies zeros for it. Its training-time provenance is not fully documented in the public release and remains an open reproducibility item.

## Native structure sources

- BioLiP3 experimental affinity source: [BioLiP affinity data](https://zhanggroup.org/BioLiP/lh_3affinity_data/00_Experiments/experiments_affinity.tsv). The pinned source snapshot hash is recorded in [NATIVE_ACCURACY_AUDIT.json](../results/NATIVE_ACCURACY_AUDIT.json). The pooled 70 includes mixed Ki/Kd endpoint types and is a post-hoc combination of two cohorts.
- OpenBind EV-A71/CVA16 2A protease: [Zenodo record 22142262](https://doi.org/10.5281/zenodo.22142262). The scored set contains 482 compounds with measured pKd from one target family; this is not broad multi-target validation.
- HiQBind source training dataset: [Figshare record 27430305](https://doi.org/10.6084/m9.figshare.27430305). The source-trained native folds use this training source; raw training labels and private records are not included here.

## Data access and redistribution

This repository does not redistribute raw affinity datasets, numerical assay labels, protected test outcomes, PDB/SDF structures or private source records. Users must obtain data from the original providers and comply with their terms. Check recorded SHA-256 values before treating a downloaded file as the pinned source snapshot.
