# Predicted-structure preparation in GeoFuse-CPI V2.0

**Status:** Implemented as an optional local research/inference preparation route. It does not change the frozen V2.0 checkpoint, the existing benchmark heads, or the V2.0 headline results.

## What it does

Given a protein sequence and a compound SMILES, the route can:

1. Obtain a sequence-matched protein structure from a user-supplied PDB/mmCIF, or from AlphaFold DB when a UniProt accession is supplied.
2. If no structure is available, optionally submit the sequence to the public ESMFold API (explicit opt-in required).
3. Generate ligand conformers from SMILES using RDKit ETKDG and MMFF (or UFF if MMFF parameters are unavailable).
4. Predict candidate binding pockets with P2Rank.
5. Prepare the protein and ligand for docking with Meeko.
6. Dock into the top-ranked pockets using AutoDock Vina.
7. Reconstruct the docked ligand and generate a 704-dimensional geometry descriptor using the V2.0 geometry feature schema.
8. Return the selected pose, docking score, pocket metadata, descriptor path, and—when model assets are available—the usual V2.0 2D task-head estimate.

The route checks protein-sequence identity/coverage against the chosen structure chain before proceeding. For predicted structures, it reports mean CA pLDDT normalized to a 0–100 scale.

## What it does not do

- It does **not** create an experimental structure. AlphaFold DB and ESMFold structures are computational predictions.
- It does **not** establish that the selected pocket binds the supplied compound.
- A P2Rank pocket probability is a pocket-ligandability score, not a compound-specific binding probability.
- An AutoDock Vina score is a docking score in kcal/mol, not pKd, Ki, Kd, or an experimental affinity.
- RDKit conformer generation is not a prediction of the bound ligand conformation.
- The generated docked-pose 704D descriptor is **not** fed into the current ATOM3D native residual expert for a 3D affinity correction. That expert was trained on experimentally bound complexes, failed its deployment gate against the pocket-only control, and has not been validated on predicted/docked poses.
- Therefore the route produces a predicted complex and a 2D model estimate, not a validated 2D+3D affinity estimate.

## Install

Use the V2.0 Python 3.11 environment:

```powershell
python -m pip install -r inference/requirements-inference.txt
python -m pip install -r inference/requirements-structure.txt
powershell -ExecutionPolicy Bypass -File scripts/install_structure_tools.ps1
python scripts/prepare_inference_assets.py
python scripts/download_esm2_weights.py
```

The installer downloads the pinned official AutoDock Vina 1.2.7 Windows executable and P2Rank 2.5.1 distribution into `%USERPROFILE%\.cache\geofuse-v2`, verifies their SHA-256 hashes, and keeps the binaries out of Git:

- Vina 1.2.7 Windows executable: `E0C4B2715E0C1A74F6E92D0F3BE0328AC97542EAFBC111E6B1EFAD897A73CCE5`
- P2Rank 2.5.1 archive: `D243F2D9036AC053FEFB9407B5FE1C85F4FE077C519FD975AC585E995FEAB274`

P2Rank requires Java 17 or later; the pipeline detects a JDK installation or respects `JAVA_HOME`.

Meeko prepares receptor and ligand PDBQT files. P2Rank predicts pockets, and AutoDock Vina performs the docking search. See the [Meeko ligand preparation guide](https://meeko.readthedocs.io/en/develop/lig_prep_basic.html), [Meeko receptor preparation guide](https://meeko.readthedocs.io/en/develop/cli_rec_prep.html), [P2Rank user guide](https://github.com/rdk/p2rank/blob/develop/documentation/user-guide.md), and [AutoDock Vina installation guide](https://autodock-vina.readthedocs.io/en/stable/installation.html).

## Run with a supplied PDB/mmCIF

```powershell
python inference/structure_pipeline.py `
  --smiles "CC(=O)Oc1ccccc1C(=O)O" `
  --fasta path/to/protein.fasta `
  --protein-pdb path/to/protein.pdb `
  --assay davis `
  --out-dir runs/aspirin_target_01
```

A matching UniProt accession can be used instead of a local structure:

```powershell
python inference/structure_pipeline.py `
  --smiles "CC(=O)Oc1ccccc1C(=O)O" `
  --fasta path/to/protein.fasta `
  --uniprot YOUR_UNIPROT_ACCESSION `
  --assay davis `
  --out-dir runs/aspirin_target_01
```

The AlphaFold DB model is used only if its chain matches the supplied sequence (minimum 95% identity and 90% query coverage). Otherwise the route reports the mismatch.

## Run when no experimental structure or UniProt accession is available

The ESMFold hosted API receives the protein sequence. This is a public external service and is rate-limited; do not use it for bulk folding or confidential sequences. The CLI requires explicit consent:

```powershell
python inference/structure_pipeline.py `
  --smiles "CC(=O)Oc1ccccc1C(=O)O" `
  --fasta path/to/protein.fasta `
  --assay davis `
  --allow-public-esmfold `
  --out-dir runs/aspirin_target_01
```

The public ESMFold route has a conservative default limit of 400 residues. Longer sequences require a matching structure or a separately provisioned local structure-prediction service. The [official ESM repository](https://github.com/facebookresearch/esm) documents the hosted fold endpoint and warns that the shared API is rate-limited and unsuitable for bulk folding.

If public sequence submission is not approved, omit `--allow-public-esmfold`. The 2D V2.0 prediction can still run; the summary will report that 3D preparation was skipped.

## Outputs

The output directory contains:

- `geofuse_structure_pipeline_summary.json`: provenance, 2D result, sequence match, pLDDT, pocket predictions, Vina scores, descriptor path, and claim boundaries.
- `source_protein.pdb` and `protein_for_docking.pdb`: downloaded/provided and protein-only structure files.
- `ligand_lowest_energy_conformer.sdf`: selected RDKit conformer.
- `receptor.pdbqt` and `ligand.pdbqt`: prepared docking inputs.
- `p2rank/`: pocket predictions.
- `docking/`: Vina poses, selected SDF pose, and per-pose 704D descriptor arrays.

The output folder is not overwritten if it already exists. Generated structures and descriptors are local research artifacts and are not committed to the public repository.

## Smoke-test record

The implementation was exercised on a public 142-residue human hemoglobin-alpha sequence with aspirin as a **software integration test only**:

| Structure source | Mean pLDDT | Top P2Rank probability | Vina score | 704D descriptor |
|---|---:|---:|---:|---|
| AlphaFold DB model | 98.06 | 0.964 | -5.945 kcal/mol | finite, shape 704 |
| ESMFold hosted API | 94.69 | 0.964 | -5.997 kcal/mol | finite, shape 704 |

These numbers only show that the structure, pocket, docking and feature-extraction stages execute. They do **not** show that aspirin binds hemoglobin, validate docking accuracy, or establish a 3D affinity improvement.

The 2D route was also exercised with the AAK1 FASTA and aspirin using the local ATOM3D-extended checkpoint; the DAVIS head returned 5.4762650089 pKd. This is a software parity check, not a measured affinity or accuracy result.

## Model-domain warning

The DAVIS/KIBA/METZ task heads have kinase-heavy benchmark source domains. The ATOM3D head is specific to the ATOM3D LBA task. None should be described as a universal pKd estimator for arbitrary proteins. The pipeline writes this warning into the JSON summary.

## External references

- [ESM repository and ESMFold API](https://github.com/facebookresearch/esm)
- [AlphaFold Protein Structure Database API](https://alphafold.ebi.ac.uk/)
- [P2Rank](https://github.com/rdk/p2rank)
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina)
- [Meeko](https://github.com/forlilab/Meeko)
