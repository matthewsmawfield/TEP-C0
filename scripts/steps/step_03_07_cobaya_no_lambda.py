#!/usr/bin/env python3
"""Step 03_07: Strict TEP Minimization without Dark Energy.

This script runs Cobaya's minimizer (PyBOBYQA/SciPy) on the joint Planck+Pantheon+
likelihoods, explicitly forcing Omega_Lambda = 0.0 in CLASS. 
This is the ultimate test of whether TEP can reconstruct the CMB without Dark Energy.
"""

import os
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "external/class/build/lib.macosx-11.1-arm64-cpython-313"))
sys.path.insert(0, str(PROJECT_ROOT))
os.environ["TEP_CLASS_PYTHONPATH"] = str(PROJECT_ROOT / "external/class/build/lib.macosx-11.1-arm64-cpython-313")

from c0_common import TEPLogger, set_step_logger, print_status, ensure_dirs, write_json, step_json_path
from cobaya.run import run

STEP_ID = "step_03_07_cobaya_no_lambda"

def run_minimization():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    ensure_dirs()
    print_status(f"Starting {STEP_ID} (Strict No-Lambda Minimization)", "TITLE")

    config = {
        "likelihood": {
            "planck_2018_highl_plik.TTTEEE": {},
            "planck_2018_lowl.TT": {},
            "planck_2018_lowl.EE": {},
            "core.pantheon_cobaya_likelihood.PantheonCobaya": {}
        },
        "theory": {
            "classy": {
                "extra_args": {
                    "N_ur": 2.0328,
                    "N_ncdm": 1,
                    "m_ncdm": 0.06,
                    "output": "tCl,pCl,lCl",
                    "l_max_scalars": 2500,
                    "lensing": "yes",
                    "tep_mode": "yes",
                    "Omega_Lambda": 0.0  # THE STRICT CONSTRAINT
                }
            }
        },
        "params": {
            "tep_epsilon_T": {"prior": {"min": -0.5, "max": 2.0}, "ref": 0.89, "proposal": 0.1, "latex": r"\epsilon_T"},
            "tep_z_T": {"prior": {"min": 0.5, "max": 15.0}, "ref": 5.0, "proposal": 0.5, "latex": r"z_T"},
            "tep_n_T": {"value": 1.0, "latex": r"n_T"},
            "H0": {"prior": {"min": 30.0, "max": 80.0}, "ref": 67.5, "proposal": 1.0, "latex": r"H_0"},
            "omega_b": {"prior": {"min": 0.018, "max": 0.030}, "ref": 0.0224, "proposal": 0.0001, "latex": r"\omega_b"},
            "omega_cdm": {"prior": {"min": 0.05, "max": 0.4}, "ref": 0.120, "proposal": 0.01, "latex": r"\omega_{cdm}"},
            "tau_reio": {"prior": {"min": 0.01, "max": 0.10}, "ref": 0.054, "proposal": 0.005, "latex": r"\tau"},
            "logA": {"prior": {"min": 2.5, "max": 3.5}, "ref": 3.044, "proposal": 0.01, "drop": True},
            "A_s": {"value": "lambda logA: 1e-10 * np.exp(logA)", "latex": r"A_s"},
            "n_s": {"prior": {"min": 0.90, "max": 1.05}, "ref": 0.966, "proposal": 0.005, "latex": r"n_s"}
        },
        "sampler": {
            "minimize": {
                "ignore_prior": False,
                "max_evals": 3000,
                "best_of": 1,
            }
        },
        "output": f"results/outputs/{STEP_ID}"
    }

    try:
        updated_info, sampler = run(config)
        mle = sampler.products()["minimum"]
        print_status(f"Minimization Complete. MLE: {mle}", "SUCCESS")
        
        # Save results
        payload = {
            "step": STEP_ID,
            "status": "completed",
            "mle": mle.to_dict() if hasattr(mle, "to_dict") else str(mle),
        }
        write_json(step_json_path(STEP_ID), payload)
        
    except Exception as e:
        print_status(f"Minimization Failed: {e}", "WARNING")
        payload = {
            "step": STEP_ID,
            "status": "failed",
            "error": str(e),
        }
        write_json(step_json_path(STEP_ID), payload)

if __name__ == "__main__":
    run_minimization()
