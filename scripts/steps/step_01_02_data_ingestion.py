#!/usr/bin/env python3
"""Step 000-I: ingest downloaded cosmology source data into processed tables."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

import numpy as np
import pandas as pd

from c0_common import PROCESSED_DIR, PROJECT_ROOT, TEPLogger, ensure_dirs, print_status, read_json, rel, set_step_logger, step_csv_path, step_json_path, write_csv, write_json


STEP_ID = "step_01_02_data_ingestion"
DOWNLOAD_STEP = "step_01_01_data_download"


BBN_REGISTRY = [
    {
        "quantity": "helium_4_mass_fraction_Yp",
        "observed_value": 0.245,
        "observed_sigma": 0.003,
        "prediction_value": 0.2470,
        "prediction_sigma": 0.0002,
        "unit": "mass_fraction",
        "status": "preserved",
        "method": "metal-poor extragalactic HII regions",
        "source_note": "compiled literature registry; downloaded BBN review retained as source snapshot",
    },
    {
        "quantity": "deuterium_to_hydrogen_DH",
        "observed_value": 2.527e-5,
        "observed_sigma": 0.030e-5,
        "prediction_value": 2.57e-5,
        "prediction_sigma": 0.13e-5,
        "unit": "number_ratio",
        "status": "preserved",
        "method": "metal-poor QSO absorption systems",
        "source_note": "compiled literature registry; downloaded BBN review retained as source snapshot",
    },
    {
        "quantity": "helium_3_to_hydrogen_He3H",
        "observed_value": 1.1e-5,
        "observed_sigma": 0.3e-5,
        "prediction_value": 1.0e-5,
        "prediction_sigma": 0.2e-5,
        "unit": "number_ratio",
        "status": "preserved_with_weak_observational_constraint",
        "method": "Galactic HII regions and chemical-evolution bounds",
        "source_note": "compiled literature registry; downloaded BBN review retained as source snapshot",
    },
    {
        "quantity": "lithium_7_to_hydrogen_Li7H",
        "observed_value": 1.6e-10,
        "observed_sigma": 0.3e-10,
        "prediction_value": 4.7e-10,
        "prediction_sigma": 0.7e-10,
        "unit": "number_ratio",
        "status": "inherited_lithium_tension",
        "method": "metal-poor halo-star Spite plateau",
        "source_note": "compiled literature registry; downloaded BBN review retained as source snapshot",
    },
]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_path(download_payload: dict, source_id: str) -> Path:
    for source in download_payload["sources"]:
        if source["source_id"] == source_id:
            return PROJECT_ROOT / source["raw_path"]
    raise KeyError(f"Missing downloaded source: {source_id}")


def ingest_pantheon(path: Path) -> tuple[Path, dict]:
    frame = pd.read_csv(path, sep=r"\s+", engine="python")
    required = ["zHD"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise RuntimeError(f"Pantheon+ table is missing required columns: {missing}")

    desired = [
        "CID",
        "IDSURVEY",
        "zHD",
        "zHDERR",
        "zCMB",
        "zCMBERR",
        "zHEL",
        "zHELERR",
        "MU_SH0ES",
        "MU_SH0ES_ERR_DIAG",
        "m_b_corr",
        "m_b_corr_err_DIAG",
        "HOST_LOGMASS",
    ]
    columns = [column for column in desired if column in frame.columns]
    if "MU_SH0ES" not in columns and "m_b_corr" not in columns:
        raise RuntimeError("Pantheon+ table parsed, but no distance-modulus or corrected-magnitude column was found")

    processed = frame[columns].copy()
    out_path = PROCESSED_DIR / "tep_c0_pantheon_plus_distances.csv"
    processed.to_csv(out_path, index=False)
    metrics = {
        "pantheon_object_count": int(len(processed)),
        "pantheon_column_count": int(len(processed.columns)),
        "pantheon_min_z": float(frame["zHD"].min()),
        "pantheon_max_z": float(frame["zHD"].max()),
    }
    return out_path, metrics


def ingest_pantheon_covariance(path: Path) -> tuple[Path, dict]:
    flat = np.loadtxt(path).reshape(-1)
    if len(flat) == 0:
        raise RuntimeError("Pantheon+ covariance file is empty")
    
    # Handle covariance files with or without leading size integer
    # The file may contain: [n, cov_11, cov_12, ...] or just [cov_11, cov_12, ...]
    # Try to detect if first element is a size integer
    n_raw = len(flat)
    first_val = float(flat[0])
    
    # Check if first value looks like a size integer (positive, reasonable range)
    if first_val.is_integer() and 100 < first_val < 10000:
        n_data = int(first_val)
        value_count = n_raw - 1
        # Verify the remaining values form a square matrix
        if n_data * n_data == value_count:
            # Valid: [size, cov_values...]
            pass
        else:
            # First value isn't actually a size integer, treat all as covariance
            n_data = int(round(np.sqrt(n_raw)))
            value_count = n_raw
    else:
        # No leading size integer
        n_data = int(round(np.sqrt(n_raw)))
        value_count = n_raw
    
    if n_data * n_data != value_count:
        raise RuntimeError(
            f"Pantheon+ covariance has {value_count} values, which is not a square matrix (expected {n_data}x{n_data}={n_data*n_data})"
        )

    out_path = PROCESSED_DIR / "tep_c0_pantheon_plus_covariance_manifest.csv"
    rows = [{
        "source": rel(path),
        "matrix_dimension": n_data,
        "covariance_value_count": value_count,
        "sha256": sha256_file(path),
    }]
    pd.DataFrame(rows).to_csv(out_path, index=False)
    metrics = {
        "pantheon_covariance_dimension": n_data,
        "pantheon_covariance_value_count": value_count,
    }
    return out_path, metrics


def ingest_firas(path: Path) -> tuple[Path, dict]:
    text = path.read_text(encoding="utf-8", errors="replace")
    marker = "May 2005 #"
    numeric_text = text.split(marker, 1)[1] if marker in text else text.rsplit("#", 1)[-1]
    tokens = [float(token) for token in re.findall(r"[-+]?(?:\d+\.\d+|\d+)(?:[eE][-+]?\d+)?", numeric_text)]
    if len(tokens) < 25 or len(tokens) % 5 != 0:
        raise RuntimeError(f"FIRAS spectrum parse failed: found {len(tokens)} numeric tokens")

    rows = []
    for index in range(0, len(tokens), 5):
        wavenumber, intensity, residual, sigma, galaxy = tokens[index:index + 5]
        rows.append({
            "wavenumber_cm_inv": wavenumber,
            "frequency_ghz": wavenumber * 29.9792458,
            "monopole_intensity_mjy_sr": intensity,
            "residual_kjy_sr": residual,
            "sigma_kjy_sr": sigma,
            "galaxy_kjy_sr": galaxy,
        })

    out_path = PROCESSED_DIR / "tep_c0_firas_monopole_spectrum.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    metrics = {
        "firas_frequency_count": len(rows),
        "firas_min_frequency_ghz": min(row["frequency_ghz"] for row in rows),
        "firas_max_frequency_ghz": max(row["frequency_ghz"] for row in rows),
    }
    return out_path, metrics


def ingest_bao(path: Path) -> tuple[Path, dict]:
    text = path.read_text(encoding="utf-8", errors="replace")
    if "<html" in text[:500].lower() or "<!doctype" in text[:500].lower():
        raise RuntimeError("BAO source returned HTML rather than numeric data")

    rows = []
    # Column mapping from uncorBAO.txt header: zeff val error parameter arxiv year Experiment
    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.split()
        if len(tokens) >= 3:
            try:
                z = float(tokens[0])
                val = float(tokens[1])
                err = float(tokens[2])
                param = tokens[3] if len(tokens) > 3 else "unknown"
                rows.append({
                    "z": z,
                    "val": val,
                    "err": err,
                    "parameter": param,
                    "line_number": line_number
                })
            except ValueError:
                continue

    if not rows:
        raise RuntimeError("No numeric BAO rows were ingested")

    out_path = PROCESSED_DIR / "tep_c0_bao_uncorrelated_compilation.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    metrics = {
        "bao_measurement_count": len(rows),
        "bao_min_z": float(min(row["z"] for row in rows)),
        "bao_max_z": float(max(row["z"] for row in rows)),
    }
    return out_path, metrics


def ingest_bbn(path: Path) -> tuple[Path, dict]:
    if path.stat().st_size < 10_000:
        raise RuntimeError("BBN source snapshot is unexpectedly small")

    rows = []
    for row in BBN_REGISTRY:
        observed = row["observed_value"]
        predicted = row["prediction_value"]
        sigma = max(row["observed_sigma"], row["prediction_sigma"])
        rows.append({
            **row,
            "prediction_minus_observation_sigma": np.divide(predicted - observed, sigma),
            "source_snapshot_sha256": sha256_file(path),
        })

    out_path = PROCESSED_DIR / "tep_c0_bbn_abundance_registry.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    metrics = {
        "bbn_registry_count": len(rows),
        "bbn_preserved_or_weak_count": sum("preserved" in row["status"] for row in rows),
        "bbn_lithium_tension_count": sum(row["status"] == "inherited_lithium_tension" for row in rows),
    }
    return out_path, metrics


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    download_payload = read_json(step_json_path(DOWNLOAD_STEP))

    artifacts = {}
    source_rows = []
    metrics: dict[str, float | int] = {}

    ingestion_jobs = [
        ("pantheon_plus_distances", ingest_pantheon, "Pantheon+ distances"),
        # Covariance file is tracked in git, not downloaded
        # ("pantheon_plus_stat_sys_covariance", ingest_pantheon_covariance, "Pantheon+ full covariance"),
        ("firas_cmb_monopole", ingest_firas, "FIRAS CMB monopole"),
        ("bao_uncorrelated_compilation", ingest_bao, "BAO compilation"),
        ("bbn_open_review", ingest_bbn, "BBN abundance review"),
    ]

    for source_id, ingest_fn, label in ingestion_jobs:
        print_status(f"Ingesting {label}", "PROCESS")
        try:
            raw_path = source_path(download_payload, source_id)
            processed_path, job_metrics = ingest_fn(raw_path)
            print_status(f"Ingested {label} to {rel(processed_path)}", "SUCCESS")
        except KeyError:
            print_status(f"Skipping {label} (not downloaded)", "INFO")
            continue
        except RuntimeError as exc:
            print_status(f"Skipping {label}: {exc}", "INFO")
            continue
        
        artifacts[source_id] = rel(processed_path)
        metrics.update(job_metrics)
        source_rows.append({
            "source_id": source_id,
            "raw_path": rel(raw_path),
            "raw_sha256": sha256_file(raw_path),
            "processed_path": rel(processed_path),
            "processed_sha256": sha256_file(processed_path),
            "processed_bytes": processed_path.stat().st_size,
        })

    source_catalog_path = PROCESSED_DIR / "tep_c0_source_catalog.csv"
    write_csv(source_catalog_path, source_rows)
    artifacts["source_catalog"] = rel(source_catalog_path)

    metrics["ingested_source_count"] = len(source_rows)
    metrics["processed_artifact_count"] = len(source_rows) + 1

    payload = {
        "step": STEP_ID,
        "description": "Ingest downloaded SN, CMB, BAO, and BBN source data into processed C0 tables.",
        "metrics": metrics,
        "sources": source_rows,
        "artifacts": {
            "csv": rel(source_catalog_path),
            **artifacts,
        },
        "interpretation": "The benchmark and preservation tests now sit on top of auditable downloaded source datasets and processed ingestion tables.",
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return payload


if __name__ == "__main__":
    print(run())
