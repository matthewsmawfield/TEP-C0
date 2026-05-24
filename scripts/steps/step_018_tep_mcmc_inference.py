#!/usr/bin/env python3
"""Step 018: independent TEP posterior diagnostic with convergence checks.

The main evidence calculation lives in step 022. This step performs an
independent MCMC pass for the mixed TEP model using the same full Pantheon+
covariance likelihood and records Gelman-Rubin / autocorrelation diagnostics.
References: Pantheon+SH0ES public distance table and covariance release.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

from c0_common import (
    TEPLogger,
    ensure_dirs,
    print_status,
    rounded,
    set_step_logger,
    step_json_path,
    write_json,
)
from step_022_three_model_comparison import (
    ModelTEP,
    PantheonData,
    RESEARCH_GRADE_RHAT_MAX,
    run_mcmc,
)

STEP_ID = "step_018_tep_mcmc_inference"


def _finite_max(values: list[float]) -> float | None:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return None
    return float(np.max(arr))


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    data = PantheonData()
    data.load()
    model = ModelTEP(pure_shear=False)

    n_walkers = int(os.getenv("TEP_STEP018_WALKERS", os.getenv("TEP_MCMC_WALKERS", "64")))
    n_steps = int(os.getenv("TEP_STEP018_STEPS", os.getenv("TEP_MCMC_STEPS", "20000")))
    burn_in = int(os.getenv("TEP_STEP018_BURN_IN", os.getenv("TEP_MCMC_BURN_IN", "5000")))

    print_status(
        f"Running mixed-model MCMC with {n_walkers} walkers, {n_steps} steps, burn-in {burn_in}",
        "PROCESS",
    )
    mcmc = run_mcmc(model, data, n_walkers=n_walkers, n_steps=n_steps, burn_in=burn_in)
    samples = np.asarray(mcmc["samples"], dtype=float)

    results = []
    for i, name in enumerate(model.param_names):
        column = samples[:, i]
        results.append({
            "parameter": name,
            "median": float(np.median(column)),
            "std": float(np.std(column)),
            "ci_16": float(np.percentile(column, 16)),
            "ci_84": float(np.percentile(column, 84)),
        })

    max_r_hat = _finite_max(mcmc["r_hat"])
    max_tau = _finite_max(mcmc["tau"])
    gelman_rubin_checked = True
    autocorr_time_checked = True
    converged = bool(
        mcmc["converged"]
        and max_r_hat is not None
        and max_r_hat < RESEARCH_GRADE_RHAT_MAX
    )

    payload = {
        "step": STEP_ID,
        "description": "Independent mixed-TEP posterior diagnostic using full Pantheon+ covariance",
        "results": results,
        "n_samples": int(len(samples)),
        "mcmc": {
            "n_walkers": n_walkers,
            "n_steps": n_steps,
            "burn_in": burn_in,
            "r_hat": mcmc["r_hat"],
            "tau": mcmc["tau"],
            "max_r_hat": max_r_hat,
            "max_tau": max_tau,
            "acceptance_fraction": mcmc["acceptance_fraction"],
            "cpu_workers": mcmc["cpu_workers"],
            "mp_context": mcmc["mp_context"],
            "autocorr_time_checked": autocorr_time_checked,
            "gelman_rubin_checked": gelman_rubin_checked,
        },
        "validation": {
            "posterior_samples_available": int(len(samples)) > 0,
            "research_grade_mcmc": converged,
            "claim_gate": "open" if converged else "blocked",
            "blockers": [] if converged else [
                "Increase MCMC length until max R-hat < 1.05 and chain length exceeds 50 autocorrelation times.",
            ],
        },
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"MCMC diagnostic completed; converged={converged}", "SUCCESS" if converged else "WARNING")
    return payload


if __name__ == "__main__":
    run()
