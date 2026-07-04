#!/usr/bin/env python3
"""Step 016: Global Likelihood Synthesis - combine all likelihoods.

References: consumes cited Pantheon+, FIRAS/Planck, BAO, and BBN step artifacts.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
from c0_common import (
    TEPLogger, ensure_dirs, print_status, read_json, 
    rounded, set_step_logger, step_json_path, write_json
)

STEP_ID = "step_03_07_likelihood_synthesis"

def run():
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Load real results - use converged Step 033 Cobaya results
    try:
        step022 = read_json(step_json_path("step_03_01_three_model_comparison"))
        step033 = read_json(step_json_path("step_03_04_cobaya_mcmc"))
    except FileNotFoundError as e:
        print_status(f"Missing dependency: {e}", "ERROR")
        raise

    m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
    logl_lcdm = step022.get('models', {}).get('M0a_LCDM', {}).get('log_likelihood_mle')
    logl_tep = step022.get('models', {}).get(m1_key, {}).get('log_likelihood_mle')
    if logl_lcdm is None or logl_tep is None:
        raise RuntimeError(f"Step 022 must report log_likelihood_mle for M0a_LCDM and {m1_key}")
    
    delta_logl = logl_tep - logl_lcdm
    
    # Use converged Step 033 parameters. Prefer the final Cobaya stats artifact
    # because the manuscript reports mean +/- std, then fall back to the chain.
    # Check both old SNe-only and new joint SNe+CMB outputs.
    cobaya_stats_paths = [
        Path("results/outputs/tep_cobaya_joint_converged_stats.txt"),  # joint SNe+CMB
        Path("results/outputs/tep_cobaya_sne_stats.txt"),             # legacy SNe-only
    ]
    cobaya_chain_paths = [
        Path("results/outputs/tep_cobaya_joint_combined.1.txt"),      # joint SNe+CMB
        Path("results/outputs/tep_cobaya_sne.1.txt"),                 # legacy SNe-only
    ]
    
    cobaya_stats_path = next((p for p in cobaya_stats_paths if p.exists()), None)
    cobaya_chain_path = next((p for p in cobaya_chain_paths if p.exists()), None)
    
    if cobaya_stats_path:
        stats = {}
        for line in cobaya_stats_path.read_text(encoding="utf-8").splitlines():
            if ":" not in line or "+/-" not in line:
                continue
            name, values = line.split(":", 1)
            mean, err = values.split("+/-", 1)
            stats[name.strip()] = (float(mean.strip()), float(err.strip()))
        h0_mcmc = stats.get("H0", (70.0, 0.0))[0]
        epsilon_t = stats.get("tep_epsilon_T", (0.1, 0.0))[0]
        mcmc_converged = step033.get("status") == "completed"
    elif cobaya_chain_path:
        chain_data = np.loadtxt(cobaya_chain_path, skiprows=1)
        # Columns: weight, minuslogpost, tep_epsilon_T, tep_z_T, H0, omega_b, omega_cdm, tau_reio, A_s, n_s, ...
        h0_samples = chain_data[:, 4]  # H0 column
        epsilon_t_samples = chain_data[:, 2]  # tep_epsilon_T column
        h0_mcmc = np.median(h0_samples)
        epsilon_t = np.median(epsilon_t_samples)
        mcmc_converged = step033.get("status") == "completed"
    else:
        # Fallback to step022 values; dimensionless models have no fitted H0
        h0_mcmc = 70.0
        epsilon_t = step022.get('models', {}).get(m1_key, {}).get('parameters_mle', {}).get('epsilon_T', 0.1)
        mcmc_converged = False
    
    # Use Planck 2018 as independent reference for tension diagnostic
    # This is a diagnostic metric, not a likelihood constraint
    h0_cmb = 67.4  # Planck 2018 fiducial (independent observational constraint)
    h0_cmb_safe = max(h0_cmb, np.finfo(float).tiny)
    trf = 1.0 / (1.0 + abs(h0_mcmc - h0_cmb)/0.5)  # Normalized agreement

    results = {
        'step': STEP_ID,
        'metrics': {
            'logl_lcdm': rounded(logl_lcdm, 1),
            'logl_tep': rounded(logl_tep, 1),
            'delta_logl': rounded(delta_logl, 1),
            'mcmc_h0': rounded(h0_mcmc, 2),
            'mcmc_ft': rounded(epsilon_t, 3),
            'tension_resolution_factor': rounded(trf, 4),
            'mcmc_converged': mcmc_converged,
        },
        'validation': {
            'claim_gate': 'open' if mcmc_converged else 'blocked',
            'blockers': [] if mcmc_converged else ['Step 033 MCMC not converged'],
            'step_033_converged': mcmc_converged,
            'cobaya_chain_available': cobaya_chain_path.exists() if 'cobaya_chain_path' in dir() else False
        },
        'interpretation': (
            f"Current Pantheon+ MLE diagnostic gives delta log L = {delta_logl:.2f} for M1 relative to M0. "
            "This is not a replacement claim unless the posterior, CMB, BBN, and robustness gates are open."
        )
    }
    
    write_json(step_json_path(STEP_ID), results)
    print_status(f"Step {STEP_ID} completed successfully", "SUCCESS")
    return results

if __name__ == "__main__":
    run()
