#!/usr/bin/env python3
import numpy as np
import shutil
import zipfile
import hashlib
import os
from pathlib import Path
import urllib.request

verbose = True
rtol = 1e-6
atol = 1e-8

## Setup
this_dir = Path(__file__).parent

zipfile_path = this_dir / "test-2D-production.zip"
url = "https://zenodo.org/records/18696440/files/test-2D-production-2026-02-19.zip"
expected_hash = "4472d27031a6fcedd30e11360814a66c075f8c14442df745e4e9d80a0dbd87f3"
expected_filenames = [
    "BOUT.restart.0.nc",
    "BOUT.restart.1.nc",
    "BOUT.restart.2.nc",
    "BOUT.restart.3.nc",
    "BOUT.restart.4.nc",
    "BOUT.restart.5.nc",
    "BOUT.restart.6.nc",
    "BOUT.restart.7.nc",
    "BOUT.restart.8.nc",
    "BOUT.restart.9.nc",
    "grid_test2_allpump.nc",
]

if verbose:
    print("2D-Production test: begin file download")

## Download files
tmp_path = zipfile_path.with_name(zipfile_path.name + ".tmp")

with urllib.request.urlopen(url, timeout=60) as response:
    if response.status != 200:
        raise RuntimeError(
            f"2D-Production test: download failed - HTTP {response.status}"
        )

    # Copy bits of the file from response to a temp file
    # This ensures no partial files are left if the download fails
    with open(tmp_path, "wb") as out_file:
        shutil.copyfileobj(response, out_file)


if verbose:
    print("2D-Production test: downloaded files")

# Rename temp file with the correct name
tmp_path.replace(zipfile_path)

with zipfile.ZipFile(zipfile_path, "r") as zf:
    zip_contents = set(zf.namelist())
    try:
        # Extract only expected grids
        for filename in expected_filenames:
            if filename in zip_contents:
                zf.extract(filename, path=this_dir)

    except Exception as e:
        print("2D-Production test: extracting test grids failed:", e)

# Check hash
with open(zipfile_path, "rb") as f:
    file_hash = hashlib.sha256(f.read()).hexdigest()
# print(file_hash)

if file_hash != expected_hash:
    raise RuntimeError(
        "2D-Production test: downloaded zip file hash does not match expected value"
    )

if verbose:
    print("2D-Production test: downloaded and extracted files")