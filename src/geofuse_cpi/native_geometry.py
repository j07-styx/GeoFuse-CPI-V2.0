"""Rigid-motion-invariant native ligand-pocket geometry features (704D)."""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

RADIAL_CENTERS = np.array([1.6, 2.2, 2.8, 3.4, 4.0, 4.7, 5.4, 6.2, 7.0, 7.8], dtype=np.float32)
ANGULAR_CENTERS = np.array([0.15, 0.5, 0.85], dtype=np.float32)
FEATURE_DIM = 49 + 49 * 10 + 49 * 3 + 10 + 3 + 3 + 2
assert FEATURE_DIM == 704


def element_class(atomic_number):
    z = np.asarray(atomic_number, dtype=np.int16)
    out = np.full(z.shape, 6, dtype=np.int16)
    out[z == 6] = 0
    out[z == 7] = 1
    out[z == 8] = 2
    out[z == 16] = 3
    out[z == 15] = 4
    out[np.isin(z, [9, 17, 35, 53])] = 5
    return out


def pair_features(graph, mode: str) -> np.ndarray:
    """Return one 704D geometry vector for a deposited bound complex.

    Modes: real native distances/angles, joint shuffled distance/angle
    assignment, or pocket-only constant-metric control. This is not docking.
    """
    if mode not in ("real", "shuffled", "pocket"):
        raise ValueError(f"Unknown geometry mode: {mode}")
    ligand_xyz = np.asarray(graph["ligand_xyz"], dtype=np.float32)
    pocket_xyz = np.asarray(graph["pocket_xyz"], dtype=np.float32)
    ligand_elements = element_class(graph["ligand_z"])
    pocket_elements = element_class(graph["pocket_z"])
    bonds = np.asarray(graph["ligand_bonds"], dtype=np.int32)
    if ligand_xyz.ndim != 2 or ligand_xyz.shape[1] != 3:
        raise ValueError("ligand_xyz must be (N,3)")
    if pocket_xyz.ndim != 2 or pocket_xyz.shape[1] != 3:
        raise ValueError("pocket_xyz must be (M,3)")
    if not len(ligand_xyz) or not len(pocket_xyz):
        raise ValueError("Empty ligand or pocket")
    if bonds.ndim != 2 or bonds.shape[0] != 2:
        raise ValueError("ligand_bonds must be (2,B)")

    members = cKDTree(pocket_xyz).query_ball_point(ligand_xyz, 8.0)
    ligand_index = np.repeat(np.arange(len(members), dtype=np.int32),
                             np.asarray([len(row) for row in members], dtype=np.int32))
    pocket_index = np.fromiter((j for row in members for j in row),
                               dtype=np.int32, count=sum(map(len, members)))
    if not len(ligand_index):
        raise ValueError("No ligand-pocket heavy-atom contacts within 8 A")
    displacement = pocket_xyz[pocket_index] - ligand_xyz[ligand_index]
    distance = np.linalg.norm(displacement, axis=1)
    if not np.isfinite(distance).all() or distance.max() > 8.001:
        raise ValueError("Invalid native contact distances")

    pair_type = ligand_elements[ligand_index] * 7 + pocket_elements[pocket_index]
    neighbors = [[] for _ in range(len(ligand_xyz))]
    for a, b in bonds.T:
        a, b = int(a), int(b)
        if not (0 <= a < len(ligand_xyz) and 0 <= b < len(ligand_xyz)) or a == b:
            continue
        neighbors[a].append(b)
        neighbors[b].append(a)

    angle = np.empty(len(distance), dtype=np.float32)
    for atom_index in np.unique(ligand_index):
        rows = np.flatnonzero(ligand_index == atom_index)
        bonded = neighbors[int(atom_index)]
        if not bonded:
            angle[rows] = 0.5
            continue
        directions = ligand_xyz[np.asarray(bonded, dtype=np.int32)] - ligand_xyz[atom_index]
        directions /= np.maximum(np.linalg.norm(directions, axis=1, keepdims=True), 1e-6)
        radial = displacement[rows] / np.maximum(distance[rows, None], 1e-6)
        angle[rows] = np.mean(np.square(radial @ directions.T), axis=1)
    angle = np.clip(angle, 0.0, 1.0)

    if mode == "shuffled":
        shift = len(distance) // 2
        distance_use = np.roll(distance, shift=shift)
        angle_use = np.roll(angle, shift=shift)
    elif mode == "pocket":
        distance_use = np.full_like(distance, 4.5)
        angle_use = np.full_like(angle, 0.5)
    else:
        distance_use, angle_use = distance, angle

    radial = np.exp(-0.5 * np.square(
        (distance_use[:, None] - RADIAL_CENTERS[None, :]) / 0.55))
    angular = np.exp(-0.5 * np.square(
        (angle_use[:, None] - ANGULAR_CENTERS[None, :]) / 0.18))
    frequency = np.bincount(pair_type, minlength=49).astype(np.float32) / len(distance)
    radial_hist = np.zeros((49, 10), dtype=np.float32)
    angular_hist = np.zeros((49, 3), dtype=np.float32)
    np.add.at(radial_hist, pair_type, radial)
    np.add.at(angular_hist, pair_type, angular)
    radial_hist /= len(distance)
    angular_hist /= len(distance)

    values = np.concatenate([
        frequency, radial_hist.ravel(), angular_hist.ravel(),
        radial.mean(axis=0), angular.mean(axis=0),
        np.log1p(np.asarray([len(ligand_xyz), len(pocket_xyz), len(distance)], dtype=np.float32)),
        np.asarray([distance_use.mean(), np.quantile(distance_use, 0.1)], dtype=np.float32),
    ]).astype(np.float32)
    if values.shape != (FEATURE_DIM,) or not np.isfinite(values).all():
        raise AssertionError("Native geometry feature must be finite 704D")
    return values
