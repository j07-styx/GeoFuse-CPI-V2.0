#!/usr/bin/env python3
"""Verify and extract the public frozen GeoFuse-CPI V2.0 inference assets."""
from __future__ import annotations

import hashlib
import json
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "inference" / "assets" / "geofuse-cpi-v2.0-inference-assets.zip"
DEST = ROOT / "models" / "v2.0"
EXPECTED_ARCHIVE_SHA256 = "71148521ae43187fcce4183cd7d388dc79daa04ebd6ace98bfacc08685a37017"
EXPECTED_FILES = {
    "GEOFUSE_CPI_V2.pt": "b30fd7f5054ea6352d53c2dc854ce0af9cd7d163970c5c4b4b0b7591772d97f4",
    "TRAIN_ANCHORS_NO_Y.npz": "42ddd4d3af1f00f15bf96985adaaa15b821fb67fe6a5f2fddf3fa796a4773769",
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> None:
    if not ARCHIVE.is_file():
        raise FileNotFoundError(f"Missing inference asset archive: {ARCHIVE}")
    archive_sha = sha256_bytes(ARCHIVE.read_bytes())
    if archive_sha != EXPECTED_ARCHIVE_SHA256:
        raise RuntimeError(f"Inference archive SHA-256 mismatch: {archive_sha}")

    DEST.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(ARCHIVE) as zf:
        manifest = json.loads(zf.read("manifest.json"))
        if manifest.get("files") != EXPECTED_FILES:
            raise RuntimeError("The embedded asset manifest differs from the pinned file hashes.")
        if set(zf.namelist()) != {"manifest.json", *EXPECTED_FILES}:
            raise RuntimeError("Unexpected file in inference asset archive.")
        for name, expected in EXPECTED_FILES.items():
            data = zf.read(name)
            actual = sha256_bytes(data)
            if actual != expected:
                raise RuntimeError(f"Asset SHA-256 mismatch for {name}: {actual}")
            target = DEST / name
            if target.exists() and sha256_bytes(target.read_bytes()) == expected:
                print(f"Already verified: {target}")
                continue
            temporary = target.with_suffix(target.suffix + ".partial")
            temporary.write_bytes(data)
            if sha256_bytes(temporary.read_bytes()) != expected:
                temporary.unlink(missing_ok=True)
                raise RuntimeError(f"Extracted asset failed verification: {name}")
            temporary.replace(target)
            print(f"Extracted and verified: {target}")

    print("PASS: frozen V2.0 checkpoint and label-free training anchors are ready.")
    print(f"Model directory: {DEST}")


if __name__ == "__main__":
    main()
