from pathlib import Path

import numpy as np
import torch
from torch import nn
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.error")


class Unified(nn.Module):
    """Frozen V2.0 ESM2/Morgan shared trunk with dataset-specific affinity heads."""

    def __init__(self):
        super().__init__()
        self.p = nn.Sequential(
            nn.Linear(1280, 128), nn.LayerNorm(128), nn.GELU(),
            nn.Dropout(.15), nn.Linear(128, 48), nn.LayerNorm(48)
        )
        self.c = nn.Sequential(
            nn.Linear(2048, 128), nn.LayerNorm(128), nn.GELU(),
            nn.Dropout(.15), nn.Linear(128, 48), nn.LayerNorm(48)
        )
        self.g = nn.Sequential(
            nn.Linear(128, 48, bias=False), nn.GELU(),
            nn.Linear(48, 48, bias=False)
        )
        self.shared = nn.Sequential(
            nn.Linear(192, 96), nn.LayerNorm(96), nn.GELU(),
            nn.Dropout(.15), nn.Linear(96, 48), nn.GELU()
        )
        self.heads = nn.ModuleDict({
            key: nn.Linear(48, 1)
            for key in ("davis", "kiba", "metz", "gensparc", "atom3d")
        })

    def representation(self, protein, compound, geometry):
        protein = self.p(protein) + self.g(geometry)
        compound = self.c(compound)
        fused = torch.cat(
            [protein, compound, protein * compound, torch.abs(protein - compound)],
            dim=1,
        )
        return self.shared(fused)

    def forward(self, protein, compound, geometry, assay):
        return self.heads[assay](self.representation(protein, compound, geometry)).squeeze(1)


def unit(values):
    return values / np.maximum(np.linalg.norm(values, axis=1, keepdims=True), 1e-8)


def fingerprints(smiles):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=True
    )
    output = []
    for source in smiles:
        molecule = Chem.MolFromSmiles(source)
        if molecule is None:
            # Match the documented ATOM3D fallback for parseable source SMILES
            # that fail strict valence sanitization (notably hypervalent P).
            molecule = Chem.MolFromSmiles(source, sanitize=False)
            if molecule is None:
                raise ValueError("SMILES cannot be parsed, even without sanitization")
            molecule.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(molecule)
        else:
            canonical = Chem.MolToSmiles(molecule, isomericSmiles=True)
            molecule = Chem.MolFromSmiles(canonical)
        fp = generator.GetFingerprint(molecule)
        array = np.zeros(2048, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, array)
        output.append(array)
    return np.asarray(output, dtype=np.float32)


class Predictor:
    def __init__(self, checkpoint, anchors):
        self.ck = torch.load(checkpoint, map_location="cpu", weights_only=True)
        self.a = np.load(anchors, allow_pickle=False)
        self.net = Unified()
        missing, unexpected = self.net.load_state_dict(self.ck["state_dict"], strict=False)
        allowed_missing = {"heads.atom3d.weight", "heads.atom3d.bias"}
        if unexpected or set(missing) - allowed_missing:
            raise ValueError(
                f"Checkpoint/model mismatch; missing={missing}, unexpected={unexpected}"
            )
        self.net.eval()
        self.available_assays = set(self.ck.get("assay_normalization", {}))

    def predict(self, assay, smiles, embeddings):
        if assay not in ("davis", "kiba", "metz", "atom3d"):
            raise ValueError("Use an assay-specific head: davis, kiba, metz, atom3d")
        if assay not in self.available_assays:
            raise ValueError(
                f"This checkpoint does not contain a trained {assay!r} head"
            )
        protein = np.asarray(embeddings, dtype=np.float32)
        compound = fingerprints(smiles)
        if protein.shape != (len(compound), 1280) or not np.isfinite(protein).all():
            raise ValueError("Expected finite ESM2 650M embeddings, 1280 values per pair")
        mean, std = self.ck["assay_normalization"][assay]
        with torch.inference_mode():
            normalized = self.net(
                torch.tensor(unit(protein) * np.sqrt(1280), dtype=torch.float32),
                torch.tensor(compound, dtype=torch.float32),
                torch.zeros((len(compound), 128), dtype=torch.float32),
                assay,
            ).numpy()
        prediction = normalized * std + mean
        if assay == "davis":
            prediction = np.maximum(prediction, 5)
        similarity = self.ck.get("similarity_heads", {}).get(assay)
        if not similarity:
            return prediction
        weight = similarity["kernel_weight"]
        if not weight:
            return prediction
        prefix = assay + "__"
        anchor_protein = self.a[prefix + "protein"]
        anchor_compound = self.a[prefix + "chemical"]
        center = self.a[prefix + "center"]
        dot = compound @ anchor_compound.T
        chemical_similarity = dot / np.maximum(
            compound.sum(1)[:, None] + anchor_compound.sum(1)[None, :] - dot, 1e-6
        )
        protein_similarity = np.clip(
            unit(protein - center) @ unit(anchor_protein - center).T, -1, 1
        )
        gamma, chemical_weight, protein_weight, _ = similarity["config"]
        kernel = np.exp(-gamma * (1 - protein_similarity))
        calibrated = similarity["mean"] + (
            chemical_similarity * kernel
            + chemical_weight * chemical_similarity
            + protein_weight * kernel
        ).astype(np.float32) @ similarity["alpha"].numpy()
        if assay == "davis":
            calibrated = np.maximum(calibrated, 5)
        return (1 - weight) * prediction + weight * calibrated
