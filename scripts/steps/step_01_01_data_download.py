#!/usr/bin/env python3
"""Step 000: download public cosmology source data for TEP-C0."""

from __future__ import annotations

import hashlib
import time
import urllib.error
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from c0_common import RAW_DIR, TEPLogger, ensure_dirs, print_status, rel, set_step_logger, step_csv_path, step_json_path, write_csv, write_json


STEP_ID = "step_01_01_data_download"
TIMEOUT_SECONDS = 300  # Increased from 60 to 300 for large covariance file downloads


SOURCE_MANIFEST = [
    {
        "source_id": "pantheon_plus_distances",
        "domain": "supernovae",
        "url": "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat",
        "target": "pantheon_plus_shoes.dat",
        "citation": "Pantheon+SH0ES DataRelease, Pantheon+ distance table.",
        "minimum_bytes": 100_000,
    },
    {
        "source_id": "pantheon_plus_stat_sys_covariance",
        "domain": "supernovae",
        "url": "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov",
        "target": "Pantheon+SH0ES.cov",
        "citation": "Pantheon+SH0ES DataRelease, full statistical plus systematic covariance matrix.",
        "minimum_bytes": 5_000_000,
    },
    {
        "source_id": "firas_cmb_monopole",
        "domain": "cmb_blackbody",
        "url": "https://lambda.gsfc.nasa.gov/data/cobe/firas/monopole_spec/firas_monopole_spec_v1.txt",
        "target": "firas_monopole_spec_v1.txt",
        "citation": "NASA LAMBDA COBE/FIRAS CMB monopole spectrum.",
        "minimum_bytes": 1_000,
    },
    {
        "source_id": "bao_uncorrelated_compilation",
        "domain": "bao",
        "url": "https://zenodo.org/records/16285883/files/uncorBAO.txt?download=1",
        "target": "uncorBAO.txt",
        "citation": "Zhao, Wang, Alam, and Ross BAO multi-survey compilation, Zenodo 10.5281/zenodo.16285883.",
        "minimum_bytes": 500,
    },
    {
        "source_id": "bbn_open_review",
        "domain": "bbn",
        "url": "https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2020.560149/full",
        "target": "bbn_frontiers_2020_review.html",
        "citation": "Pizzone et al. 2020 Frontiers review of BBN reaction rates and abundance constraints.",
        "minimum_bytes": 10_000,
    },
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, target: Path) -> tuple[str, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "TEP-C0-data-downloader/0.1"})
    temp_target = target.with_suffix(target.suffix + ".tmp")
    
    # Stream download with progress reporting for large files
    with urllib.request.urlopen(request, timeout=TIMEOUT_SECONDS) as response:
        total_size = int(response.headers.get('Content-Length', 0))
        downloaded = 0
        chunk_size = 8192
        
        with temp_target.open('wb') as f:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                f.write(chunk)
                downloaded += len(chunk)
                
                # Report progress for large files (>1MB)
                if total_size > 1_000_000 and downloaded % (1024 * 1024) == 0:
                    progress = (downloaded / total_size * 100) if total_size > 0 else 0
                    print_status(f"Download progress: {progress:.1f}% ({downloaded/1024/1024:.1f}MB/{total_size/1024/1024:.1f}MB)", "INFO")
    
    temp_target.replace(target)
    return "downloaded", ""


def download_source(source: dict) -> dict:
    """Download a single source and return its metadata."""
    target_path = RAW_DIR / source["target"]
    status = "missing"
    error = ""

    print_status(f"Processing source: {source['source_id']}", "PROCESS")
    try:
        status, error = download(source["url"], target_path)
        print_status(f"Downloaded {source['source_id']}", "SUCCESS")
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        error = str(exc)
        if target_path.exists() and target_path.stat().st_size >= source["minimum_bytes"]:
            status = "cached"
            print_status(f"Using cached {source['source_id']}", "SUCCESS")
        else:
            print_status(f"Failed to download {source['source_id']}: {error}", "ERROR")
            return {
                "source": source,
                "status": "failed",
                "error": error,
                "exception": str(exc)
            }

    size_bytes = target_path.stat().st_size
    if size_bytes < source["minimum_bytes"]:
        print_status(f"Source {source['source_id']} too small: {size_bytes} bytes", "ERROR")
        return {
            "source": source,
            "status": "failed",
            "error": f"File too small: {size_bytes} bytes"
        }

    return {
        "source": source,
        "status": status,
        "error": error,
        "size_bytes": size_bytes,
        "target_path": target_path
    }


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    rows = []
    downloaded_count = 0
    cached_count = 0
    total_bytes = 0

    # Use parallel downloads for efficiency
    print_status(f"Downloading {len(SOURCE_MANIFEST)} sources in parallel", "INFO")
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_to_source = {executor.submit(download_source, source): source for source in SOURCE_MANIFEST}
        
        for future in as_completed(future_to_source):
            result = future.result()
            source = result["source"]
            
            if result["status"] == "failed":
                error = result.get("error", "Unknown error")
                print_status(f"Failed to download {source['source_id']}: {error}", "ERROR")
                raise RuntimeError(f"Failed to download {source['source_id']} from {source['url']}: {error}")
            
            status = result["status"]
            target_path = result["target_path"]
            size_bytes = result["size_bytes"]
            
            if status == "downloaded":
                downloaded_count += 1
            elif status == "cached":
                cached_count += 1
            total_bytes += size_bytes

            rows.append({
                "source_id": source["source_id"],
                "domain": source["domain"],
                "status": status,
                "bytes": size_bytes,
                "sha256": sha256_file(target_path),
                "raw_path": rel(target_path),
                "url": source["url"],
                "citation": source["citation"],
                "download_epoch": int(time.time()),
                "error": result.get("error", ""),
            })

    manifest_csv = step_csv_path(STEP_ID)
    write_csv(manifest_csv, rows)

    payload = {
        "step": STEP_ID,
        "description": "Download and cache public cosmology source datasets for the TEP-C0 empirical front end.",
        "metrics": {
            "source_count": len(SOURCE_MANIFEST),
            "downloaded_source_count": downloaded_count,
            "cached_source_count": cached_count,
            "total_raw_bytes": total_bytes,
        },
        "sources": rows,
        "artifacts": {
            "csv": rel(manifest_csv),
            **{row["source_id"]: row["raw_path"] for row in rows},
        },
        "interpretation": "The C0 pipeline now begins from explicit public data products rather than benchmark-only synthetic inputs.",
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
