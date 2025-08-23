#!/usr/bin/env python3
"""
Simple helper that runs the ready‑to‑use fetch script.
Assumes the Bash script is in the repo root.
"""

import os
import subprocess
from pathlib import Path
import shutil

DATA_DIR = Path("data")

def run_fetcher() -> Path:
    """Call the Bash script, return absolute path to output bundle """
    script = Path("./data_fetcher.sh").resolve()
    if not script.is_file():
        raise FileNotFoundError(f"{script} not found")
    print(f"Running data downloader script: {script}")
    subprocess.run([str(script)], check=True)
    zipfile = next(Path(".").glob("SlAREB1_sequences_bundle_SGN_first.zip"))
    if not zipfile:
        raise RuntimeError("Missing expected zip file.")
    extracted = DATA_DIR / "bundle"
    if extracted.exists():
        shutil.rmtree(extracted)
    extracted.mkdir(parents=True, exist_ok=True)
    print(f"Extracting {zipfile} -> {extracted}")
    subprocess.run(["unzip", "-q", str(zipfile), "-d", str(extracted)], check=True)
    return extracted

if __name__ == "__main__":
    run_fetcher()