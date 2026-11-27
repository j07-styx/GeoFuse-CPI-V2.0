# Frozen V2.0 inference package

## 2D raw-input prediction

The raw predictor accepts one SMILES string and one protein FASTA sequence. It computes the ESM-2 650M embedding, Morgan fingerprint and the frozen V2.0 assay-head output.

Install and prepare:

    python -m pip install torch==2.14.0
    python -m pip install -r inference/requirements-inference.txt
    python scripts/prepare_inference_assets.py
    python scripts/download_esm2_weights.py

Run:

    python inference/predict.py --smiles "CC(=O)Oc1ccccc1C(=O)O" --fasta examples/AAK1.fasta --assay davis

Assay choices are davis, kiba and metz. DAVIS and METZ scores are pKd; KIBA is its composite benchmark score. No native 3D geometry is used in this raw-input route.

## Native expert

The native residual expert is a separate model component. It consumes prepared x2d=(N,423) and real/shuffled/pocket=(N,704) arrays. It does not predict or dock poses from raw PDB/SDF inputs.

    python -m pip install -r inference/requirements-native.txt
    python inference/predict_native_prepared.py --features prepared_noY.npz --out predictions.csv

The input archive must be label-free and must not contain affinity/target arrays.

## Reference environment

- Python 3.11.9
- torch 2.14.0
- numpy 2.4.6
- RDKit 2025.9.2
- fair-esm 2.0.0
- lightgbm 4.7.0 for native expert inference

ESM-2 weights are downloaded separately and SHA-256 verified by scripts/download_esm2_weights.py.
