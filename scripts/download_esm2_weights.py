#!/usr/bin/env python3
"""Download the public ESM-2 650M weights and verify the pinned SHA-256."""
from __future__ import annotations

import argparse
import hashlib
import urllib.request
from pathlib import Path

URL = "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt"
EXPECTED_SHA256 = "ea9d0522b335a8778dea6535a65301f10208dece28cd5865482b0b1fc446168c"
EXPECTED_SIZE = 2_604_537_549
DEFAULT_OUTPUT = Path("models/esm2_t33_650M_UR50D.pt")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    if output.is_file():
        digest = sha256_file(output)
        if output.stat().st_size == EXPECTED_SIZE and digest == EXPECTED_SHA256:
            print(f"Already verified: {output}")
            return
        raise RuntimeError(f"Existing ESM file is not the pinned checkpoint: {output}")

    partial = output.with_suffix(output.suffix + ".partial")
    partial.unlink(missing_ok=True)
    request = urllib.request.Request(URL, headers={"User-Agent": "GeoFuse-CPI-V2.0-reproducibility/1.0"})
    print(f"Downloading {URL}")
    try:
        with urllib.request.urlopen(request, timeout=120) as response, partial.open("wb") as target:
            length = response.headers.get("Content-Length")
            if length and int(length) != EXPECTED_SIZE:
                raise RuntimeError(f"Unexpected ESM checkpoint size from source: {length}")
            total = 0
            while True:
                block = response.read(8 << 20)
                if not block:
                    break
                total += len(block)
                if total > EXPECTED_SIZE:
                    raise RuntimeError("ESM download exceeded the pinned size.")
                target.write(block)
                print(f"Downloaded {total:,}/{EXPECTED_SIZE:,} bytes", flush=True)
        if partial.stat().st_size != EXPECTED_SIZE:
            raise RuntimeError(f"Incomplete ESM download: {partial.stat().st_size} bytes")
        digest = sha256_file(partial)
        if digest != EXPECTED_SHA256:
            raise RuntimeError(f"ESM SHA-256 mismatch: {digest}")
        partial.replace(output)
    except Exception:
        partial.unlink(missing_ok=True)
        raise

    print(f"PASS: ESM-2 650M weights verified at {output}")


if __name__ == "__main__":
    main()
