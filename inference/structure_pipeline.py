#!/usr/bin/env python3
"""GeoFuse-CPI V2.0 predicted-structure and docking preparation pipeline.

This pipeline always keeps the V2.0 2D prediction separate from predicted
protein structures and docked ligand poses. It can generate a 704D geometry
descriptor from a predicted complex, but it deliberately does NOT apply the
current ATOM3D native-geometry residual expert: that expert was trained on
experimentally bound complexes and is not deployment-eligible.

Structure sources, in order:
  1. Explicit user-supplied PDB.
  2. AlphaFold DB, when a UniProt accession is supplied and sequence-matched.
  3. ESMFold hosted API, for supported sequence lengths.

Pocket prediction: P2Rank. Ligand conformers: RDKit ETKDG + MMFF/UFF.
Docking: AutoDock Vina via its official Windows executable. Ligand/receptor
preparation: Meeko. All generated structures are computational predictions,
not experimental observations or validated binding poses.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import requests
import gemmi
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBIO, PDBParser, Select, is_aa
from Bio.SeqUtils import seq1
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "inference"))
sys.path.insert(0, str(ROOT / "src"))

from geofuse_portable_kinase import Predictor
from geofuse_raw_kinase import Encoder, clean as clean_sequence
from geofuse_cpi.native_geometry import pair_features

try:
    from meeko import (
        MoleculePreparation,
        PDBQTWriterLegacy,
        PDBQTMolecule,
        Polymer,
        RDKitMolCreate,
    )
except ImportError as exc:
    raise SystemExit(
        "Structure inference requires Meeko. Install "
        "inference/requirements-structure.txt in the V2.0 Python environment."
    ) from exc

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
DEFAULT_P2RANK_HOME = Path.home() / ".cache" / "geofuse-v2" / "p2rank_2.5.1"
DEFAULT_VINA = Path.home() / ".cache" / "geofuse-v2" / "vina.exe"
UNITS = {"davis": "pKd", "kiba": "KIBA score", "metz": "pKd", "atom3d": "pK"}
AA_CUSTOM = {"MSE": "M", "SEC": "U", "PYL": "O"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_query_sequence(sequence: str | None, fasta: Path | None) -> str:
    if sequence:
        raw = sequence
    elif fasta:
        raw = fasta.read_text(encoding="utf-8")
    else:
        raise ValueError("Provide --sequence or --fasta")
    return clean_sequence(raw)


def chain_sequences(pdb_path: Path) -> list[dict[str, Any]]:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("query", str(pdb_path))
    model = next(structure.get_models())
    result = []
    for chain in model:
        residues = [
            residue for residue in chain.get_residues()
            if residue.id[0] in (" ", "H_MSE")
            and is_aa(residue, standard=False)
            and "CA" in residue
        ]
        seq = "".join(seq1(residue.resname, custom_map=AA_CUSTOM) or "X" for residue in residues)
        if seq:
            result.append({"chain": chain.id, "sequence": seq, "residues": residues})
    if not result:
        raise ValueError(f"No protein chains with CA atoms found in {pdb_path}")
    return result


def alignment_quality(query: str, target: str) -> tuple[float, float]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -3.0
    aligner.extend_gap_score = -0.25
    alignment = aligner.align(query, target)[0]
    query_blocks, target_blocks = alignment.aligned
    identical = 0
    aligned = 0
    query_covered = 0
    for (q0, q1), (t0, t1) in zip(query_blocks, target_blocks):
        length = min(int(q1 - q0), int(t1 - t0))
        identical += sum(query[int(q0) + i] == target[int(t0) + i] for i in range(length))
        aligned += length
        query_covered += int(q1 - q0)
    identity = identical / aligned if aligned else 0.0
    coverage = min(1.0, query_covered / max(1, len(query)))
    return float(identity), float(coverage)


def select_matching_chain(
    sequence: str, pdb_path: Path, requested_chain: str | None,
    min_identity: float, min_coverage: float,
) -> dict[str, Any]:
    candidates = chain_sequences(pdb_path)
    if requested_chain is not None:
        candidates = [row for row in candidates if row["chain"] == requested_chain]
        if not candidates:
            raise ValueError(f"Chain {requested_chain!r} is absent from {pdb_path}")
    scored = []
    for row in candidates:
        identity, coverage = alignment_quality(sequence, row["sequence"])
        scored.append({**row, "identity": identity, "coverage": coverage})
    best = max(scored, key=lambda row: (row["identity"] * row["coverage"], row["coverage"]))
    if best["identity"] < min_identity or best["coverage"] < min_coverage:
        raise ValueError(
            "Structure sequence does not adequately match the supplied protein "
            f"(chain={best['chain']!r}, identity={best['identity']:.3f}, "
            f"query_coverage={best['coverage']:.3f}; required "
            f"{min_identity:.2f}/{min_coverage:.2f})."
        )
    return best


def normalize_structure_file(source: Path, destination: Path) -> None:
    suffix = source.suffix.lower()
    if suffix in (".cif", ".mmcif"):
        structure = gemmi.read_structure(str(source))
        structure.write_pdb(str(destination))
    else:
        shutil.copy2(source, destination)
    if not destination.is_file() or destination.stat().st_size < 100:
        raise RuntimeError(f"Could not normalize structure file: {source}")


def fetch_afdb_structure(accession: str, out_path: Path, timeout: int) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_-]{4,16}", accession):
        raise ValueError("UniProt accession contains unsupported characters")
    response = requests.get(
        AFDB_API.format(accession=accession.upper()),
        timeout=timeout,
        headers={"User-Agent": "GeoFuse-CPI-V2.0-structure-preparation/1.0"},
    )
    response.raise_for_status()
    records = response.json()
    if not isinstance(records, list) or not records:
        raise RuntimeError(f"AlphaFold DB returned no model for {accession}")
    model = next((row for row in records if row.get("pdbUrl")), None)
    if model is None:
        raise RuntimeError(f"AlphaFold DB response has no PDB coordinate URL for {accession}")
    pdb_response = requests.get(model["pdbUrl"], timeout=timeout)
    pdb_response.raise_for_status()
    if "ATOM" not in pdb_response.text[:20000]:
        raise RuntimeError("AlphaFold DB returned a file without PDB ATOM records")
    out_path.write_text(pdb_response.text, encoding="utf-8")
    return str(model["pdbUrl"])


def fold_with_esmfold(sequence: str, out_path: Path, timeout: int, max_length: int) -> str:
    if len(sequence) > max_length:
        raise RuntimeError(
            f"ESMFold hosted API limit configured at {max_length} residues; "
            f"received {len(sequence)}. Supply a matching PDB/UniProt structure "
            "or use a separately provisioned local structure-prediction service."
        )
    response = requests.post(
        ESMFOLD_URL,
        data=sequence,
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "User-Agent": "GeoFuse-CPI-V2.0-structure-preparation/1.0",
        },
        timeout=timeout,
    )
    response.raise_for_status()
    text = response.text
    if "ATOM  " not in text and "ATOM " not in text:
        raise RuntimeError(
            "ESMFold endpoint did not return PDB coordinates: " + text[:240].replace("\n", " ")
        )
    out_path.write_text(text, encoding="utf-8")
    return ESMFOLD_URL


class ProteinOnlySelect(Select):
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_model(self, model):
        return int(model.id) == 0

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " " and is_aa(residue, standard=True)

    def accept_atom(self, atom):
        element = (atom.element or atom.get_name()[0]).upper()
        return element not in ("H", "D") and atom.get_altloc() in (" ", "A")


def write_protein_only(source: Path, destination: Path, chain_id: str) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(source))
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(str(destination), ProteinOnlySelect(chain_id))
    if not destination.is_file() or destination.stat().st_size < 100:
        raise RuntimeError("Failed to write a protein-only PDB for docking")


def mean_plddt(pdb_path: Path, chain_id: str) -> float | None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("confidence", str(pdb_path))
    values = []
    for residue in structure[0][chain_id].get_residues():
        if is_aa(residue, standard=False) and "CA" in residue:
            values.append(float(residue["CA"].bfactor))
    if not values:
        return None
    mean = float(np.mean(values))
    # ESMFold's hosted PDB stores pLDDT on a 0-1 scale; AlphaFold DB
    # stores it on a 0-100 scale. Normalize both to the documented 0-100.
    if max(values) <= 1.5:
        mean *= 100.0
    return float(min(100.0, max(0.0, mean)))


def acquire_protein_structure(
    sequence: str, out_dir: Path, protein_pdb: Path | None, uniprot: str | None,
    requested_chain: str | None, structure_kind: str, timeout: int,
    esmfold_max_length: int, allow_public_esmfold: bool,
) -> dict[str, Any]:
    source_path = out_dir / "source_protein.pdb"
    source = None
    source_url = None
    kind = structure_kind

    if protein_pdb is not None:
        if not protein_pdb.is_file():
            raise FileNotFoundError(protein_pdb)
        normalize_structure_file(protein_pdb, source_path)
        source = "user_supplied_structure"
    elif uniprot:
        try:
            source_url = fetch_afdb_structure(uniprot, source_path, timeout)
            select_matching_chain(sequence, source_path, requested_chain, 0.95, 0.90)
            source = "alphafold_database"
            kind = "predicted"
        except Exception as exc:
            # A UniProt isoform/construct mismatch should not silently become
            # the structure for the supplied sequence.
            source_path.unlink(missing_ok=True)
            afdb_error = str(exc)
            if not allow_public_esmfold:
                raise RuntimeError(
                    "The AlphaFold DB model did not match the supplied sequence. "
                    "A fallback would send the sequence to the public ESMFold API; "
                    "rerun with --allow-public-esmfold to permit that."
                ) from exc
            source_url = fold_with_esmfold(sequence, source_path, timeout, esmfold_max_length)
            source = "esmfold_api_fallback"
            kind = "predicted"
    else:
        if not allow_public_esmfold:
            raise RuntimeError(
                "No local/AlphaFold DB structure was supplied. Public ESMFold folding "
                "sends the protein sequence to an external service; pass "
                "--allow-public-esmfold to opt in. The 2D predictor remains available."
            )
        source_url = fold_with_esmfold(sequence, source_path, timeout, esmfold_max_length)
        source = "esmfold_api"
        kind = "predicted"

    min_identity, min_coverage = (0.95, 0.90) if kind == "predicted" else (0.85, 0.70)
    match = select_matching_chain(
        sequence, source_path, requested_chain, min_identity, min_coverage
    )
    clean_path = out_dir / "protein_for_docking.pdb"
    write_protein_only(source_path, clean_path, match["chain"])
    confidence = mean_plddt(source_path, match["chain"]) if kind == "predicted" else None

    result = {
        "source": source,
        "source_url": source_url,
        "structure_kind": kind,
        "source_file": str(source_path),
        "protein_only_file": str(clean_path),
        "selected_chain": match["chain"],
        "sequence_identity": match["identity"],
        "query_sequence_coverage": match["coverage"],
        "mean_pLDDT": confidence,
        "pLDDT_note": (
            "Mean CA B-factor interpreted as pLDDT for predicted models."
            if kind == "predicted" else
            "Not reported: experimental B-factors are not pLDDT."
        ),
    }
    if "afdb_error" in locals():
        result["alphafold_db_fallback_reason"] = afdb_error
    return result


def find_java_home() -> str | None:
    def valid_home(candidate: Path) -> bool:
        executable = candidate / "bin" / ("java.exe" if os.name == "nt" else "java")
        if not executable.is_file():
            return False
        try:
            result = subprocess.run(
                [str(executable), "-version"], capture_output=True,
                text=True, timeout=10,
            )
            version_text = result.stdout + result.stderr
            match = re.search(r'(?:version\s+"|openjdk\s+)(\d+)', version_text, re.I)
            return bool(match and int(match.group(1)) >= 17)
        except Exception:
            return False

    configured = os.environ.get("JAVA_HOME")
    if configured and valid_home(Path(configured)):
        return configured
    if os.name == "nt":
        roots = [
            Path("C:/Program Files/Java"),
            Path("C:/Program Files/Eclipse Adoptium"),
            Path("C:/Program Files/Microsoft"),
        ]
        candidates = []
        for root in roots:
            if root.is_dir():
                candidates.extend(root.glob("jdk-*"))
        candidates.sort(
            key=lambda path: tuple(int(v) for v in re.findall(r"\d+", path.name)),
            reverse=True,
        )
        for candidate in candidates:
            if valid_home(candidate):
                return str(candidate)
    java = shutil.which("java")
    if java:
        resolved = Path(java).resolve()
        if resolved.parent.name.lower() == "bin" and valid_home(resolved.parent.parent):
            return str(resolved.parent.parent)
    return None


def locate_p2rank(home_arg: Path | None) -> tuple[Path, Path]:
    home = home_arg or Path(os.environ.get("GEOFUSE_P2RANK_HOME", DEFAULT_P2RANK_HOME))
    if home.is_file():
        home = home.parent
    exe = home / ("prank.bat" if os.name == "nt" else "prank")
    if not exe.is_file():
        raise FileNotFoundError(
            f"P2Rank launcher not found at {exe}. Run scripts/install_structure_tools.ps1 "
            "or set GEOFUSE_P2RANK_HOME."
        )
    return home, exe


def run_p2rank(
    protein_pdb: Path, output_dir: Path, structure_kind: str,
    p2rank_home_arg: Path | None, timeout: int,
) -> list[dict[str, Any]]:
    home, exe = locate_p2rank(p2rank_home_arg)
    output_dir.mkdir(parents=True, exist_ok=True)
    args = [str(exe), "predict", "-f", str(protein_pdb), "-o", str(output_dir),
            "-visualizations", "0"]
    if structure_kind == "predicted":
        # P2Rank's AlphaFold profile avoids treating model B-factors as
        # experimental temperature factors.
        args.extend(["-c", "alphafold"])
    env = os.environ.copy()
    java_home = find_java_home()
    if os.name == "nt":
        if not java_home:
            raise RuntimeError("P2Rank requires Java 17+; set JAVA_HOME to a JDK installation.")
        env["JAVA_HOME"] = java_home
        command = subprocess.list2cmdline(args)
        completed = subprocess.run(
            ["cmd.exe", "/d", "/c", command],
            cwd=str(home), env=env, capture_output=True, text=True, timeout=timeout,
        )
    else:
        completed = subprocess.run(
            args, cwd=str(home), env=env, capture_output=True, text=True, timeout=timeout,
        )
    if completed.returncode != 0:
        raise RuntimeError(
            "P2Rank failed: " + (completed.stderr + "\n" + completed.stdout)[-3000:]
        )
    csv_files = sorted(output_dir.glob("*_predictions.csv"))
    if not csv_files:
        raise RuntimeError("P2Rank completed but produced no *_predictions.csv file")
    pockets = []
    with csv_files[0].open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        for raw in reader:
            row = {(key or "").strip(): (value or "").strip()
                   for key, value in raw.items()}
            try:
                pockets.append({
                    "rank": int(row["rank"]),
                    "score": float(row["score"]),
                    "probability": float(row["probability"]),
                    "center": [
                        float(row["center_x"]), float(row["center_y"]), float(row["center_z"])
                    ],
                    "residue_ids": row.get("residue_ids", "").split(),
                })
            except (KeyError, ValueError):
                continue
    pockets.sort(key=lambda row: (row["probability"], row["score"]), reverse=True)
    if not pockets:
        raise RuntimeError("P2Rank found no usable pocket predictions")
    return pockets


def generate_ligand_conformer(smiles: str, out_dir: Path, seed: int, n_conformers: int = 12):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse the supplied SMILES")
    molecule = Chem.AddHs(molecule)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useRandomCoords = False
    conformer_ids = list(AllChem.EmbedMultipleConfs(
        molecule, numConfs=int(n_conformers), params=params
    ))
    if not conformer_ids:
        raise RuntimeError("RDKit ETKDG failed to generate a ligand conformer")

    if AllChem.MMFFHasAllMoleculeParams(molecule):
        energies = AllChem.MMFFOptimizeMoleculeConfs(molecule, maxIters=500)
        method = "MMFF"
    else:
        energies = AllChem.UFFOptimizeMoleculeConfs(molecule, maxIters=500)
        method = "UFF"
    best_index = min(range(len(energies)), key=lambda index: energies[index][1])
    best_conf_id = conformer_ids[best_index]
    conformer_path = out_dir / "ligand_lowest_energy_conformer.sdf"
    writer = Chem.SDWriter(str(conformer_path))
    writer.write(molecule, confId=best_conf_id)
    writer.close()

    setup = MoleculePreparation().prepare(molecule, conformer_id=best_conf_id)[0]
    pdbqt_string, success, error_message = PDBQTWriterLegacy.write_string(
        setup, add_index_map=True
    )
    if not success or not pdbqt_string:
        raise RuntimeError("Meeko ligand preparation failed: " + error_message)
    ligand_pdbqt = out_dir / "ligand.pdbqt"
    ligand_pdbqt.write_text(pdbqt_string, encoding="utf-8")
    return {
        "molecule": molecule,
        "pdbqt": ligand_pdbqt,
        "conformer_sdf": conformer_path,
        "conformer_method": method,
        "conformer_count": len(conformer_ids),
        "selected_conformer_energy": float(energies[best_index][1]),
        "protonation_policy": "As supplied in input SMILES; no pH/tautomer enumeration.",
    }


def prepare_receptor(protein_pdb: Path, out_dir: Path) -> Path:
    polymer = Polymer.from_pdb_file(str(protein_pdb))
    polymer.parameterize(MoleculePreparation())
    rigid_pdbqt, _flex_pdbqt = PDBQTWriterLegacy.write_string_from_polymer(polymer)
    if not rigid_pdbqt:
        raise RuntimeError("Meeko generated an empty receptor PDBQT")
    receptor_path = out_dir / "receptor.pdbqt"
    receptor_path.write_text(rigid_pdbqt, encoding="utf-8")
    return receptor_path


def parse_vina_score(pose_path: Path) -> float | None:
    for line in pose_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "REMARK VINA RESULT:" in line:
            try:
                return float(line.split("REMARK VINA RESULT:", 1)[1].split()[0])
            except (ValueError, IndexError):
                return None
    return None


def pose_geometry_descriptor(
    pose_path: Path, protein_pdb: Path, output_npy: Path,
) -> tuple[np.ndarray, Path]:
    pdbqt = PDBQTMolecule.from_file(str(pose_path), name="GeoFuse_docked_pose", poses_to_read=1)
    reconstructed = RDKitMolCreate.from_pdbqt_mol(pdbqt)
    if not reconstructed:
        raise RuntimeError("Meeko could not reconstruct the docked ligand bond graph")
    ligand = Chem.RemoveHs(reconstructed[0])
    if ligand.GetNumConformers() == 0:
        raise RuntimeError("Docked ligand has no 3D coordinates")
    ligand_xyz = np.asarray(ligand.GetConformer().GetPositions(), dtype=np.float32)
    ligand_z = np.asarray([atom.GetAtomicNum() for atom in ligand.GetAtoms()], dtype=np.int16)
    edges = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in ligand.GetBonds()]
    ligand_bonds = (
        np.asarray(edges, dtype=np.int32).T if edges
        else np.empty((2, 0), dtype=np.int32)
    )

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(protein_pdb))
    protein_atoms = []
    for residue in structure[0].get_residues():
        if not is_aa(residue, standard=True):
            continue
        for atom in residue.get_atoms():
            element = (atom.element or atom.get_name()[0]).upper()
            if element in ("H", "D") or atom.get_altloc() not in (" ", "A"):
                continue
            protein_atoms.append(atom)
    if not protein_atoms:
        raise RuntimeError("No standard protein heavy atoms available for geometry extraction")
    protein_xyz = np.asarray([atom.coord for atom in protein_atoms], dtype=np.float32)
    periodic = Chem.GetPeriodicTable()
    protein_z = []
    for atom in protein_atoms:
        element = (atom.element or atom.get_name()[0]).title()
        try:
            protein_z.append(periodic.GetAtomicNumber(element))
        except Exception:
            protein_z.append(0)
    protein_z = np.asarray(protein_z, dtype=np.int16)

    tree = cKDTree(protein_xyz)
    neighborhoods = tree.query_ball_point(ligand_xyz, r=8.0)
    nearby_indices = sorted({index for row in neighborhoods for index in row})
    if len(nearby_indices) > 2000:
        nearest_distance = cKDTree(ligand_xyz).query(protein_xyz[nearby_indices])[0]
        keep = np.argsort(nearest_distance)[:2000]
        nearby_indices = [nearby_indices[int(index)] for index in keep]
    if len(nearby_indices) < 15:
        raise RuntimeError(
            f"Docked pose has only {len(nearby_indices)} receptor heavy atoms within 8 A"
        )
    graph = {
        "ligand_xyz": ligand_xyz,
        "ligand_z": ligand_z,
        "ligand_bonds": ligand_bonds,
        "pocket_xyz": protein_xyz[nearby_indices],
        "pocket_z": protein_z[nearby_indices],
    }
    feature = pair_features(graph, "real")
    if feature.shape != (704,) or not np.isfinite(feature).all():
        raise RuntimeError("Generated geometry descriptor failed the 704D finite-value check")
    np.save(output_npy, feature.astype(np.float32))
    pose_sdf = output_npy.with_name("selected_docked_ligand_pose.sdf")
    Chem.MolToMolFile(ligand, str(pose_sdf))
    return feature, pose_sdf


def dock_top_pockets(
    receptor_pdbqt: Path, ligand_pdbqt: Path, pockets: list[dict[str, Any]],
    protein_pdb: Path, output_dir: Path, vina_path: Path, pocket_count: int,
    box_size: float, exhaustiveness: int, cpu: int, timeout: int,
) -> list[dict[str, Any]]:
    if not vina_path.is_file():
        raise FileNotFoundError(
            f"AutoDock Vina executable not found at {vina_path}. "
            "Run scripts/install_structure_tools.ps1 or set GEOFUSE_VINA."
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for pocket in pockets[:pocket_count]:
        rank = pocket["rank"]
        pose_path = output_dir / f"pocket_{rank}_poses.pdbqt"
        center_x, center_y, center_z = pocket["center"]
        command = [
            str(vina_path), "--receptor", str(receptor_pdbqt),
            "--ligand", str(ligand_pdbqt),
            "--center_x", str(center_x), "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(box_size), "--size_y", str(box_size),
            "--size_z", str(box_size),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", "5", "--energy_range", "3",
            "--cpu", str(cpu), "--out", str(pose_path),
        ]
        try:
            completed = subprocess.run(
                command, capture_output=True, text=True, timeout=timeout
            )
            if completed.returncode != 0 or not pose_path.is_file():
                results.append({
                    "pocket_rank": rank, "status": "docking_failed",
                    "error": (completed.stderr + "\n" + completed.stdout)[-1800:],
                })
                continue
            score = parse_vina_score(pose_path)
            row = {
                "pocket_rank": rank,
                "p2rank_score": pocket["score"],
                "p2rank_probability": pocket["probability"],
                "pocket_center_A": pocket["center"],
                "vina_score_kcal_per_mol": score,
                "pose_pdbqt": str(pose_path),
                "status": "docked",
            }
            descriptor_path = output_dir / f"pocket_{rank}_pose_704d.npy"
            try:
                feature, pose_sdf = pose_geometry_descriptor(
                    pose_path, protein_pdb, descriptor_path
                )
                row.update({
                    "geometry_feature_dim": int(feature.shape[0]),
                    "geometry_feature_file": str(descriptor_path),
                    "pose_sdf": str(pose_sdf),
                    "geometry_status": "generated_from_predicted_docked_pose",
                })
            except Exception as exc:
                row["geometry_status"] = "feature_extraction_failed"
                row["geometry_error"] = str(exc)
            results.append(row)
        except subprocess.TimeoutExpired:
            results.append({
                "pocket_rank": rank, "status": "docking_timeout",
                "error": f"Vina exceeded {timeout} seconds",
            })
    return results


def run_2d_prediction(
    smiles: str, sequence: str, assay: str, checkpoint: Path,
    anchors: Path, esm_checkpoint: Path,
) -> dict[str, Any]:
    if not checkpoint.is_file():
        raise FileNotFoundError(f"GeoFuse checkpoint not found: {checkpoint}")
    if not anchors.is_file():
        raise FileNotFoundError(f"GeoFuse training anchors not found: {anchors}")
    if not esm_checkpoint.is_file():
        raise FileNotFoundError(
            f"ESM-2 weights not found: {esm_checkpoint}. Run scripts/download_esm2_weights.py."
        )
    predictor = Predictor(checkpoint, anchors)
    try:
        encoder = Encoder(esm_checkpoint)
        embedding = encoder.encode(sequence)
        value = float(predictor.predict(assay, [smiles], embedding[None, :])[0])
    finally:
        predictor.a.close()
    return {
        "status": "ok",
        "assay_head": assay,
        "prediction": value,
        "unit": UNITS[assay],
        "checkpoint": str(checkpoint),
        "checkpoint_sha256": sha256_file(checkpoint),
        "model_domain_warning": (
            "DAVIS/KIBA/METZ are benchmark task heads with kinase-heavy source "
            "domains; ATOM3D is specific to the ATOM3D LBA task. These outputs "
            "are not universal pKd estimates for arbitrary proteins."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    query = parser.add_mutually_exclusive_group(required=True)
    query.add_argument("--sequence", help="Protein amino-acid sequence.")
    query.add_argument("--fasta", type=Path, help="Single-protein FASTA file.")
    parser.add_argument("--smiles", required=True)
    parser.add_argument("--assay", required=True, choices=tuple(UNITS))
    parser.add_argument("--protein-pdb", type=Path)
    parser.add_argument("--uniprot")
    parser.add_argument("--chain")
    parser.add_argument("--structure-kind", choices=("experimental", "predicted"), default="experimental")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--checkpoint", type=Path, default=None)
    parser.add_argument("--anchors", type=Path, default=ROOT / "models" / "v2.0" / "TRAIN_ANCHORS_NO_Y.npz")
    parser.add_argument("--esm-checkpoint", type=Path, default=ROOT / "models" / "esm2_t33_650M_UR50D.pt")
    parser.add_argument("--p2rank-home", type=Path, default=None)
    parser.add_argument("--vina", type=Path, default=None)
    parser.add_argument("--pocket-count", type=int, default=3)
    parser.add_argument("--box-size", type=float, default=22.0)
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20261127)
    parser.add_argument("--structure-timeout", type=int, default=300)
    parser.add_argument("--esmfold-max-length", type=int, default=400)
    parser.add_argument(
        "--allow-public-esmfold", action="store_true",
        help="Explicitly permit sending the supplied protein sequence to the public ESMFold API.",
    )
    parser.add_argument("--docking-timeout", type=int, default=240)
    parser.add_argument("--skip-2d", action="store_true",
                        help="Structure-preparation smoke test only; no affinity prediction.")
    args = parser.parse_args()

    if args.pocket_count < 1 or args.box_size < 12 or args.exhaustiveness < 1:
        parser.error("pocket-count/exhaustiveness must be positive and box-size >= 12 A")
    sequence = read_query_sequence(args.sequence, args.fasta)
    if args.protein_pdb is None and not args.uniprot and args.structure_kind == "experimental":
        # The default for a sequence-only request is predicted, not experimental.
        args.structure_kind = "predicted"

    checkpoint = args.checkpoint
    if checkpoint is None:
        local_extension = ROOT / "models" / "v2.0" / "GEOFUSE_CPI_V2_ATOM3D_LOCAL.pt"
        public_frozen = ROOT / "models" / "v2.0" / "GEOFUSE_CPI_V2.pt"
        checkpoint = local_extension if local_extension.is_file() else public_frozen
    vina_path = args.vina or Path(os.environ.get("GEOFUSE_VINA", DEFAULT_VINA))

    if args.out_dir.exists():
        parser.error(f"Output directory already exists; refusing to overwrite: {args.out_dir}")
    args.out_dir.mkdir(parents=True, exist_ok=False)

    summary: dict[str, Any] = {
        "model": "GeoFuse-CPI V2.0",
        "status": "started",
        "input": {
            "sequence_length": len(sequence),
            "smiles": args.smiles,
            "assay": args.assay,
        },
        "2d_prediction": {"status": "skipped"} if args.skip_2d else {},
        "structure_pipeline": {},
        "3d_affinity_prediction": {
            "status": "not_applied",
            "prediction": None,
            "reason": (
                "The current ATOM3D 704D residual expert was trained on experimentally "
                "bound complexes, failed its deployment gate against the pocket-only "
                "control, and has not been validated on predicted/docked poses. The "
                "generated pose and descriptor are therefore structural research outputs, "
                "not a calibrated 3D affinity prediction."
            ),
        },
        "notices": [
            "The public ESMFold API receives the protein sequence only when --allow-public-esmfold is explicitly supplied.",
            "ESMFold/AlphaFold coordinates are computationally predicted structures, not experimental structures.",
            "P2Rank probability estimates pocket ligandability, not whether this compound binds.",
            "AutoDock Vina score is a docking score in kcal/mol, not pKd/Ki/Kd or an experimental affinity.",
            "RDKit conformer generation does not determine the bound ligand conformation.",
            "No experimental 3D superiority or docking accuracy is claimed by this inference route.",
        ],
    }

    if not args.skip_2d:
        try:
            summary["2d_prediction"] = run_2d_prediction(
                args.smiles, sequence, args.assay, checkpoint,
                args.anchors, args.esm_checkpoint,
            )
        except Exception as exc:
            summary["2d_prediction"] = {"status": "failed", "error": str(exc)}

    try:
        structure = acquire_protein_structure(
            sequence, args.out_dir, args.protein_pdb, args.uniprot,
            args.chain, args.structure_kind, args.structure_timeout,
            args.esmfold_max_length, args.allow_public_esmfold,
        )
        pockets = run_p2rank(
            Path(structure["protein_only_file"]),
            args.out_dir / "p2rank",
            structure["structure_kind"],
            args.p2rank_home,
            args.structure_timeout,
        )
        ligand = generate_ligand_conformer(
            args.smiles, args.out_dir, args.seed
        )
        receptor_pdbqt = prepare_receptor(
            Path(structure["protein_only_file"]), args.out_dir
        )
        docked = dock_top_pockets(
            receptor_pdbqt, ligand["pdbqt"], pockets,
            Path(structure["protein_only_file"]),
            args.out_dir / "docking", vina_path,
            args.pocket_count, args.box_size, args.exhaustiveness,
            args.cpu, args.docking_timeout,
        )
        successful = [
            row for row in docked
            if row.get("status") == "docked" and row.get("vina_score_kcal_per_mol") is not None
        ]
        best = min(successful, key=lambda row: row["vina_score_kcal_per_mol"]) if successful else None
        summary["structure_pipeline"] = {
            "status": "complete" if best else "partial_no_successful_docking",
            "protein_structure": structure,
            "ligand_conformer": {
                key: (str(value) if isinstance(value, Path) else value)
                for key, value in ligand.items() if key != "molecule"
            },
            "pocket_method": "P2Rank 2.5.1",
            "pockets": pockets,
            "docking_method": "AutoDock Vina 1.2.7",
            "docking_runs": docked,
            "selected_docking_result": best,
            "selected_704d_geometry_file": best.get("geometry_feature_file") if best else None,
            "vina_executable": str(vina_path),
            "p2rank_home": str(locate_p2rank(args.p2rank_home)[0]),
        }
    except Exception as exc:
        summary["structure_pipeline"] = {
            "status": "failed",
            "error": str(exc),
        }

    two_d_ok = summary["2d_prediction"].get("status") == "ok"
    structure_ok = summary["structure_pipeline"].get("status") == "complete"
    if args.skip_2d and structure_ok:
        summary["status"] = "structure_preparation_only"
    elif two_d_ok and structure_ok:
        summary["status"] = "2d_prediction_and_predicted_complex_prepared"
    elif two_d_ok:
        summary["status"] = "2d_prediction_available_structure_pipeline_incomplete"
    elif structure_ok:
        summary["status"] = "predicted_complex_prepared_2d_prediction_unavailable"
    else:
        summary["status"] = "partial_or_failed"

    summary_path = args.out_dir / "geofuse_structure_pipeline_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, allow_nan=False))
    print(f"SUMMARY_FILE={summary_path}")
    return 0 if two_d_ok or structure_ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
