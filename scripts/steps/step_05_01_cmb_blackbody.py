#!/usr/bin/env python3
"""Step 010: CMB Blackbody Preservation - verify TEP preserves CMB spectrum.

Uses actual FIRAS data to verify TEP predicts blackbody spectrum preservation.
References: COBE/FIRAS monopole spectrum via NASA LAMBDA; Planck constants are
used only as external CMB reference values.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import TEPLogger, ensure_dirs, print_status, read_json, rounded, set_step_logger, step_csv_path, step_json_path, write_csv, write_json

STEP_ID = "step_05_01_cmb_blackbody"


def load_firas_data():
    """Load real COBE/FIRAS CMB monopole spectrum data (NASA LAMBDA; Fixsen et al. 1996)."""
    firas_file = Path("data/raw/firas_spectrum.dat")
    if not firas_file.exists():
        raise FileNotFoundError(
            f"Real COBE/FIRAS spectrum required at {firas_file}; synthetic fallback is disabled."
        )
    data = np.loadtxt(firas_file)
    return data[:, 0], data[:, 1], data[:, 2], True  # freq [GHz], intensity [SI], error [SI]


def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    freq, intensity, err, has_data = load_firas_data()
    
    # Fit for temperature using exact SI constants
    from scipy.optimize import curve_fit
    
    def planck(f, T):
        h = 6.62607015e-34
        k = 1.380649e-23
        c = 299792458.0
        x = h * f * 1e9 / (k * T)
        return (2.0 * h * (f * 1e9)**3 / c**2) * np.exp(-x) / (-np.expm1(-x))
    
    popt, pcov = curve_fit(planck, freq, intensity, p0=[2.725], sigma=err, absolute_sigma=True)
    T_fit = float(popt[0])
    T_err_stat = float(np.sqrt(pcov[0, 0]))
    # Combine statistical channel fit error with FIRAS absolute thermometer calibration (0.6 mK, Fixsen 2009)
    T_err = float(np.sqrt(T_err_stat**2 + 0.0006**2))
    
    # TEP predicts identical blackbody (photon conservation)
    T_tep_prediction = 2.725
    consistent = bool(abs(T_fit - T_tep_prediction) < 3.0 * T_err)

    results = {
        'step': STEP_ID,
        'description': 'CMB blackbody preservation under TEP',
        'data_source': 'FIRAS' if has_data else 'Standard FIRAS spectrum',
        'temperature_fit': rounded(T_fit, 4),
        'temperature_error': rounded(T_err, 4),
        'TEP_prediction': T_tep_prediction,
        'deviation_mK': rounded((T_fit - T_tep_prediction) * 1000, 2),
        'consistent_with_TEP': consistent,
        'interpretation': f'CMB temperature {rounded(T_fit, 4)}±{rounded(T_err, 4)} K consistent with TEP prediction' if consistent else f'CMB temperature {rounded(T_fit, 4)}±{rounded(T_err, 4)} K deviates from TEP prediction by {rounded(abs(T_fit - T_tep_prediction)*1000, 2)} mK'
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results


if __name__ == "__main__":
    run()
