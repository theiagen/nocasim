#!/usr/bin/env python3
"""Download HuNoV probe sequences from Bhamidipati 2025 Supplementary File 2."""

import sys
import urllib.request
from pathlib import Path

URL = (
    "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12216758/"
    "bin/41598_2025_14659_MOESM2_ESM.txt"
)


def main():
    outdir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/probes")
    outdir.mkdir(parents=True, exist_ok=True)
    dest = outdir / "hunov_probes.fasta"

    print(f"Downloading probes from {URL}")
    urllib.request.urlretrieve(URL, dest)
    print(f"Saved to {dest}")


if __name__ == "__main__":
    main()
