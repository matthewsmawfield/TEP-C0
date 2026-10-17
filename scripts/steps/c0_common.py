#!/usr/bin/env python3
"""Shared utilities for the TEP-C0 benchmark pipeline."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable
import numpy as np

try:
    from scripts.utils.logger import TEPLogger, print_status, set_step_logger
except ImportError:
    class TEPLogger:
        def __init__(self, *args, **kwargs): pass
    def print_status(msg, status="INFO"): print(f"[{status}] {msg}")
    def set_step_logger(logger): pass

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"

def load_model_comparison() -> dict[str, Any]:
    """Verified current SN comparison in the existing reporting schema."""
    path = RESULTS_DIR / 'current_model_comparison.json'
    if not path.exists():
        raise FileNotFoundError('Run scripts/steps/refresh_observer_comparison.py before reporting SN evidence')
    current = json.loads(path.read_text())
    if current.get('redshift') != 'zHD' or current.get('observer_redshift') != 'zHEL':
        raise ValueError('Current SN comparison has an unverified redshift convention')
    for rel, expected in current['source_sha256'].items():
        if hashlib.sha256((PROJECT_ROOT/rel).read_bytes()).hexdigest() != expected:
            raise ValueError(f'Comparison source changed: {rel}; regenerate current comparison')
    if hashlib.sha256((RAW_DIR/'Pantheon+SH0ES.cov').read_bytes()).hexdigest() != current['covariance_sha256']:
        raise ValueError('Comparison covariance changed; regenerate current comparison')
    aliases = {'M1_zT5': 'M1_NoLambda_zT5', 'M1_zT100': 'M1_Unscreened_zT100'}
    models, bayes = {}, {}
    baseline = current['models']['M0a_LCDM']['logZ']
    for name, record in current['models'].items():
        assert np.isfinite([record['chi2'], record['log_likelihood'], record['logZ'], record['logZ_error']]).all()
        assert record['logZ_error'] >= 0
        assert np.isclose(record['BIC'], record['k']*np.log(current['sample_size'])-2*record['log_likelihood'], atol=1e-8, rtol=0)
        assert np.isclose(record['BF_vs_LCDM'], np.exp(record['logZ']-baseline), rtol=1e-10, atol=0)
        key = aliases.get(name, name)
        models[key] = dict(record, parameters_mle=record['parameters'],
                           log_evidence=record['logZ'], log_evidence_error=record['logZ_error'],
                           log_likelihood_mle=record['log_likelihood'])
        bayes[f'BF_{key}_vs_M0a_LCDM'] = record['BF_vs_LCDM']
        bayes[f'ln_BF_{key}_vs_M0a_LCDM'] = record['logZ']-baseline
    return {'source_file': str(path), 'status': 'completed',
            'data': {'n_supernovae': current['sample_size']}, 'models': models,
            'bayes_factors': bayes, 'conventions': {k:v for k,v in current.items() if k not in ['models', 'source_sha256']},
            'validation': {'research_grade_provenance': True, 'evidence_numerically_verified': True,
                           'model_comparison_gate': 'open', 'blockers': [],
                           'scope': 'Matched SN evidence calculation; not full cosmological observational closure.'}}


C_KM_S = 299_792.458

# Cosmological constants for transport kernel and distance calculations
# Used by step_02_01_transport_kernel and related steps
H0_KM_S_MPC = 70.0
OMEGA_M = 0.3
OMEGA_L = 0.7
OMEGA_K = 0.0

def ensure_dirs() -> None:
    for d in [RAW_DIR, RESULTS_DIR, FIGURES_DIR, PROCESSED_DIR]:
        d.mkdir(parents=True, exist_ok=True)

def e_z(z: np.ndarray | float) -> np.ndarray | float:
    zp1 = 1 + np.asarray(z)
    # Add radiation term for completeness (important at high z)
    OMEGA_R = 9.0e-5  # Radiation density parameter (photons + neutrinos)
    return np.sqrt(OMEGA_M * zp1 ** 3 + OMEGA_R * zp1 ** 4 + OMEGA_K * zp1 ** 2 + OMEGA_L)

def h_km_s_mpc(z: np.ndarray | float) -> np.ndarray | float:
    return H0_KM_S_MPC * e_z(z)

def comoving_distance_mpc(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    h_values = np.maximum(h_km_s_mpc(z), np.finfo(float).tiny)
    integrand = np.divide(C_KM_S, h_values)
    distances = np.zeros_like(z)
    if len(z) > 1:
        dz = np.diff(z)
        trapezoids = 0.5 * (integrand[1:] + integrand[:-1]) * dz
        distances[1:] = np.cumsum(trapezoids)
    return distances

def angular_diameter_distance_mpc(z: np.ndarray) -> np.ndarray:
    return comoving_distance_mpc(z) / (1 + np.asarray(z))

class TEPEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer, np.floating, np.bool_)):
            return obj.item()
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, cls=TEPEncoder) + "\n", encoding="utf-8")

def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))

def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows: return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

def rel(path: Path) -> str:
    return str(path.relative_to(PROJECT_ROOT))

def step_json_path(step_id: str) -> Path:
    return RESULTS_DIR / f"{step_id}.json"

def step_csv_path(step_id: str) -> Path:
    return RESULTS_DIR / f"{step_id}.csv"

def figure_path(step_id: str) -> Path:
    return FIGURES_DIR / f"{step_id}.png"

def rounded(value: float, digits: int = 6) -> float:
    return round(float(value), digits)
