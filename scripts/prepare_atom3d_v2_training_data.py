#!/usr/bin/env python3
"""Prepare ATOM3D LBA official TRAIN only for the GeoFuse-CPI V2.0 head.

The output is a local research cache, not for redistribution. The official
validation/test directories are never opened. No label-based filtering occurs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from geofuse_cpi.native_geometry import pair_features
from atom3d_lba_reader import TrainOnlyLMDB

RDLogger.DisableLog("rdApp.error")
AA3 = {
    "ALA":"A", "ARG":"R", "ASN":"N", "ASP":"D", "CYS":"C",
    "GLU":"E", "GLN":"Q", "GLY":"G", "HIS":"H", "ILE":"I",
    "LEU":"L", "LYS":"K", "MET":"M", "PHE":"F", "PRO":"P",
    "SER":"S", "THR":"T", "TRP":"W", "TYR":"Y", "VAL":"V",
}
FPGEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048, includeChirality=True
)
PERIODIC = Chem.GetPeriodicTable()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def protein_sequence(item) -> str:
    protein = item.get("atoms_protein")
    if not isinstance(protein, dict) or not all(k in protein for k in ("columns", "data")):
        raise ValueError("Original protein atom table unavailable")
    col = {name: i for i, name in enumerate(protein["columns"])}
    required = ("name", "resname", "chain", "residue", "insertion_code")
    if any(name not in col for name in required):
        raise ValueError("Original protein atom table lacks CA sequence fields")
    chains = defaultdict(list)
    seen = set()
    for record in protein["data"]:
        if record[col["name"]] != "CA":
            continue
        aa = AA3.get(str(record[col["resname"]]).upper())
        if aa is None:
            continue
        chain = str(record[col["chain"]])
        key = (chain, str(record[col["residue"]]), str(record[col["insertion_code"]]))
        if key in seen:
            continue
        seen.add(key)
        chains[chain].append(aa)
    sequences = ["".join(values) for values in chains.values() if len(values) >= 30]
    if not sequences:
        raise ValueError("Could not reconstruct an observed CA protein chain")
    sequence = max(sequences, key=len)
    if len(sequence) > 2048:
        raise ValueError(f"Observed protein chain exceeds V2.0 ESM input limit: {len(sequence)}")
    return sequence


def morgan(smiles: str) -> tuple[np.ndarray, str]:
    mol = Chem.MolFromSmiles(smiles)
    mode = "rdkit_sanitized"
    if mol is None:
        # Some deposited source SMILES are parseable but fail strict valence
        # sanitization (e.g. hypervalent P). Preserve the source graph and use
        # a documented unsanitized Morgan fallback; do not invent chemistry.
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            raise ValueError("Source SMILES cannot be parsed even without sanitization")
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(mol)
        mode = "rdkit_unsanitized_valence_fallback"
    fp = FPGEN.GetFingerprint(mol)
    out = np.zeros(2048, dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, out)
    if out.shape != (2048,):
        raise AssertionError("Morgan fingerprint shape mismatch")
    return out.astype(np.float32), mode


def geometry(item) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ligand = item["atoms_ligand"]
    pocket = item["atoms_pocket"]
    bonds = item["bonds"]
    ligand_elements = ligand["element"].astype(str).str.upper().to_numpy()
    pocket_elements = pocket["element"].astype(str).str.upper().to_numpy()
    ligand_keep = np.flatnonzero(~np.isin(ligand_elements, ("H", "D")))
    pocket_keep = np.flatnonzero(~np.isin(pocket_elements, ("H", "D")))
    if not (5 <= len(ligand_keep) <= 200 and 15 <= len(pocket_keep) <= 2000):
        raise ValueError("Heavy-atom ligand/pocket size outside fixed feature protocol")
    ligand_xyz = ligand[["x", "y", "z"]].to_numpy(dtype=np.float32)[ligand_keep]
    pocket_xyz = pocket[["x", "y", "z"]].to_numpy(dtype=np.float32)[pocket_keep]
    def atomic_number(symbol):
        normalized = symbol[:1].upper() + symbol[1:].lower()
        try:
            return PERIODIC.GetAtomicNumber(normalized)
        except Exception:
            return 0
    ligand_z = np.asarray([atomic_number(x) for x in ligand_elements[ligand_keep]], dtype=np.int16)
    pocket_z = np.asarray([atomic_number(x) for x in pocket_elements[pocket_keep]], dtype=np.int16)
    old_to_new = np.full(len(ligand), -1, dtype=np.int32)
    old_to_new[ligand_keep] = np.arange(len(ligand_keep), dtype=np.int32)
    edges = []
    for a, b in bonds[["atom1", "atom2"]].to_numpy(dtype=np.int32):
        if 0 <= a < len(old_to_new) and 0 <= b < len(old_to_new):
            aa, bb = int(old_to_new[a]), int(old_to_new[b])
            if aa >= 0 and bb >= 0 and aa != bb:
                edges.append((aa, bb))
    bond_index = np.asarray(edges, dtype=np.int32).T if edges else np.empty((2, 0), dtype=np.int32)
    graph = {
        "ligand_xyz": ligand_xyz,
        "ligand_z": ligand_z,
        "ligand_bonds": bond_index,
        "pocket_xyz": pocket_xyz,
        "pocket_z": pocket_z,
    }
    return tuple(pair_features(graph, mode) for mode in ("real", "shuffled", "pocket"))


def group_ids(sequences: list[str]) -> np.ndarray:
    from sklearn.feature_extraction.text import CountVectorizer
    x = CountVectorizer(analyzer="char", ngram_range=(3, 3),
                        binary=True, dtype=np.int32).fit_transform(sequences)
    x.data[:] = 1
    sizes = np.asarray(x.sum(axis=1)).ravel()
    overlap = (x @ x.T).tocoo()
    similarity = overlap.data.astype(np.float64) / np.maximum(
        1.0, sizes[overlap.row] + sizes[overlap.col] - overlap.data
    )
    keep = (overlap.row < overlap.col) & (similarity >= 0.28)
    rows = np.r_[np.arange(len(sequences), dtype=np.int32), overlap.row[keep], overlap.col[keep]]
    cols = np.r_[np.arange(len(sequences), dtype=np.int32), overlap.col[keep], overlap.row[keep]]
    graph = csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                       shape=(len(sequences), len(sequences)))
    return connected_components(graph, directed=False)[1].astype(np.int32)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-lmdb", required=True, type=Path)
    parser.add_argument("--folds-json", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        parser.error(f"Refusing to overwrite existing prepared data: {args.out}")
    if not args.folds_json.is_file():
        parser.error("Frozen TRAIN-only four-fold assignment JSON is missing")

    folds_doc = json.loads(args.folds_json.read_text(encoding="utf-8"))
    if folds_doc.get("status") != "TRAIN_ONLY_HEURISTIC_PROTEIN_GROUPED_4FOLD":
        parser.error("Unexpected or non-TRAIN-only fold assignment")
    if folds_doc.get("n_train") != 3507 or folds_doc.get("original_VAL_or_TEST_accessed") is not False:
        parser.error("Fold assignment does not match the official TRAIN-only contract")
    fold_id = np.full(3507, -1, dtype=np.int8)
    for fold, indices in enumerate(folds_doc["folds"]):
        for index in indices:
            if not 0 <= int(index) < 3507 or fold_id[int(index)] != -1:
                parser.error("Fold indices are invalid or repeated")
            fold_id[int(index)] = fold
    if np.any(fold_id < 0) or [int(np.sum(fold_id == i)) for i in range(4)] != [877, 877, 877, 876]:
        parser.error("Frozen fold coverage/sizes changed")

    db = TrainOnlyLMDB(args.train_lmdb)
    n = len(db)
    ids, smiles, sequences, labels, fps, modes = [], [], [], [], [], []
    real, shuffled, pocket = [], [], []
    try:
        for i in range(n):
            item = db.get(i)
            pdb_id = str(item["id"])
            smi = str(item["smiles"])
            y = float(item["scores"]["neglog_aff"])
            if not np.isfinite(y):
                raise ValueError(f"Nonfinite official TRAIN label at source row {i}")
            seq = protein_sequence(item)
            fp, mode = morgan(smi)
            r, s, p = geometry(item)
            ids.append(pdb_id)
            smiles.append(smi)
            sequences.append(seq)
            labels.append(y)
            fps.append(fp)
            modes.append(mode)
            real.append(r)
            shuffled.append(s)
            pocket.append(p)
            if (i + 1) % 250 == 0 or i + 1 == n:
                print(f"ATOM3D V2.0 TRAIN preparation {i+1}/{n}", flush=True)
    finally:
        db.close()

    group = group_ids(sequences)
    if len(np.unique(group)) != 1004 or int(np.bincount(group).max()) != 302:
        raise AssertionError("Protein 3-mer grouping no longer matches frozen source audit")
    for g in np.unique(group):
        if len(np.unique(fold_id[group == g])) != 1:
            raise AssertionError("A heuristic protein group crosses the frozen four folds")

    data = {
        "source_index": np.arange(n, dtype=np.int32),
        "ids": np.asarray(ids),
        "smiles": np.asarray(smiles),
        "sequences": np.asarray(sequences),
        "y": np.asarray(labels, dtype=np.float32),
        "morgan": np.stack(fps).astype(np.float32),
        "fingerprint_mode": np.asarray(modes),
        "real": np.stack(real).astype(np.float32),
        "shuffled": np.stack(shuffled).astype(np.float32),
        "pocket": np.stack(pocket).astype(np.float32),
        "fold_id": fold_id,
        "protein_group_id": group,
    }
    if data["morgan"].shape != (3507, 2048) or any(data[k].shape != (3507, 704) for k in ("real", "shuffled", "pocket")):
        raise AssertionError("Prepared ATOM3D V2.0 feature shapes are incorrect")
    if not all(np.isfinite(data[k]).all() for k in ("y", "morgan", "real", "shuffled", "pocket")):
        raise AssertionError("Prepared ATOM3D V2.0 features contain nonfinite values")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("xb") as stream:
        np.savez_compressed(stream, **data)
    report = {
        "status": "ATOM3D_V2_0_TRAIN_ONLY_FEATURES_PREPARED",
        "n_official_train": n,
        "official_val_or_test_opened": False,
        "official_test_rescored": False,
        "source_train_lmdb_sha256": sha256(args.train_lmdb / "data.mdb"),
        "fold_assignment_sha256": sha256(args.folds_json),
        "prepared_npz_sha256": sha256(args.out),
        "protein_sequence_count_unique": len(set(sequences)),
        "protein_group_count_3mer_jaccard_028": int(len(np.unique(group))),
        "fold_sizes": [int(np.sum(fold_id == i)) for i in range(4)],
        "morgan_sanitized_count": int(np.sum(np.asarray(modes) == "rdkit_sanitized")),
        "morgan_unsanitized_valence_fallback_count": int(np.sum(np.asarray(modes) == "rdkit_unsanitized_valence_fallback")),
        "morgan_unparsed_count": 0,
        "geometry_feature_dim": 704,
        "license": "CC BY-NC-ND; local noncommercial research only; do not redistribute dataset derivatives or trained weights without permission",
        "limits": [
            "Official 466-row validation and 490-row test are not used in this integration run.",
            "Four-fold grouping is based on observed CA-chain 3-mer Jaccard >=0.28, not aligned sequence identity.",
            "Ligand stereochemistry is only as represented by the deposited source SMILES; unsanitized fallback is explicitly flagged.",
            "Native geometry features require the experimentally bound ligand-pocket complex and are not pose generation or docking."
        ],
    }
    report_path = args.out.with_suffix(".audit.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
