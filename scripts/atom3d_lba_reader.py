"""Read-only decoder for the official ATOM3D LBA TRAIN LMDB only.

This helper refuses non-TRAIN paths and requires exactly 3,507 records.
It never opens the official validation or test directories.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import lmdb
import msgpack
import pandas as pd


class TrainOnlyLMDB:
    def __init__(self, path: str | Path):
        self.path = Path(path).resolve()
        normalized = self.path.as_posix().lower()
        if self.path.name.lower() != "train" or "/official_raw_train/" not in normalized:
            raise ValueError("Only the explicit official_raw_train/.../data/train LMDB is allowed")
        if not self.path.exists():
            raise FileNotFoundError(self.path)
        self.env = lmdb.open(str(self.path), readonly=True, lock=False,
                             max_readers=2, readahead=False, meminit=False)
        with self.env.begin() as txn:
            self.size = int(txn.get(b"num_examples"))
            self.format = txn.get(b"serialization_format").decode("utf-8")
        if self.size != 3507:
            self.close()
            raise ValueError(f"Expected 3,507 official TRAIN records, got {self.size}")
        if self.format not in ("json", "msgpack"):
            self.close()
            raise ValueError(f"Refusing unsafe/unrecognized LMDB format: {self.format}")

    def __len__(self):
        return self.size

    def get(self, index: int):
        if not 0 <= index < self.size:
            raise IndexError(index)
        with self.env.begin() as txn:
            raw = txn.get(str(index).encode("ascii"))
        if raw is None:
            raise ValueError(f"Missing official TRAIN record {index}")
        payload = gzip.decompress(raw)
        item = json.loads(payload) if self.format == "json" else msgpack.unpackb(payload, raw=False)
        for key in ("atoms_pocket", "atoms_ligand", "bonds"):
            value = item.get(key)
            if isinstance(value, dict) and all(k in value for k in ("columns", "index", "data")):
                item[key] = pd.DataFrame(**value)
        return item

    def close(self):
        if getattr(self, "env", None) is not None:
            self.env.close()
            self.env = None
