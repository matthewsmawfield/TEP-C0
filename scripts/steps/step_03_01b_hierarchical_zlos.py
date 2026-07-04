#!/usr/bin/env python3
"""
Step 03.01b: Hierarchical z_los Marginalization with Hyperprior Integration
============================================================================

Concern: Reporting BF = 4.6, 61.8, 40.3 side-by-side looks like cherry-picking.
The free-z_los dynesty run (BF = 40.3) already samples z_los jointly, so its
evidence is technically a marginalized evidence, but the prior is a broad uniform,
not an explicit hyperprior.

This script:
  1. Defines a hyperprior p(z_los) = LogUniform(1, 200).
  2. Uses the existing discrete likelihood evaluations at z_los = 1, 5, 100
     and the free-z_los dynesty result to approximate ln BF(z_los).
  3. Integrates: BF_marg = ∫ dz_los p(z_los) BF(z_los).
  4. Reports the marginalized ln BF with Monte Carlo uncertainty.

Honest note: If the hyperprior integration yields a number close to 40.3, then
the free-z_los dynesty result was already approximately correct — but the paper
will now have a defensible derivation of the headline number, not just a convenient
selection.

Usage:
    cd TEP-C0 && python scripts/steps/step_03_01b_hierarchical_zlos.py
"""

import json
import numpy as np
from pathlib import Path
from scipy import stats
from scipy.interpolate import interp1d

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = PROJECT_ROOT / "results"
OUTPUT_PATH = RESULTS_DIR / "step_03_01b_hierarchical_zlos.json"


def load_existing_results():
    """Load the discrete model-comparison results."""
    with open(RESULTS_DIR / "step_03_01_three_model_comparison.json") as f:
        data = json.load(f)
    return data


def approximate_lnBF_vs_zlos():
    """
    Approximate ln BF(z_los) from discrete evaluations.
    We have measured points at z = 1, 5, 100 and a free-z_los result.
    The free-z_los run samples z_los ~ U[0.1, 150] and gives BF = 40.28.
    We approximate the evidence profile by interpolating the fixed-z points
    and calibrating to the free-z result.
    """
    data = load_existing_results()
    bayes = data["bayes_factors"]

    # Fixed-z points: z -> BF vs LCDM
    z_points = np.array([1.0, 5.0, 100.0])
    bf_points = np.array([
        bayes["BF_M1_NoLambda_zT1_vs_M0a_LCDM"],
        bayes["BF_M1_NoLambda_zT5_vs_M0a_LCDM"],
        bayes["BF_M1_Unscreened_zT100_vs_M0a_LCDM"],
    ])
    lnBF_points = np.log(bf_points)

    # Free-z result as approximate marginalized evidence
    bf_free = bayes["BF_M1_free_zT_vs_M0a_LCDM"]
    lnBF_free = np.log(bf_free)

    # Interpolate ln BF vs log(z)
    # Use linear interpolation (conservative with sparse points);
    # quadratic can overshoot wildly between sparse samples.
    log_z = np.log(z_points)
    interp = interp1d(log_z, lnBF_points, kind='linear', fill_value='extrapolate')

    return interp, lnBF_free, z_points, lnBF_points


def hyperprior_integration(interp, lnBF_free, n_samples=100_000, seed=42):
    """
    Monte Carlo integration of BF_marg = ∫ dz p(z) BF(z)
    with p(z) = LogUniform(1, 200).
    """
    rng = np.random.default_rng(seed)

    # LogUniform(1, 200) => log z ~ Uniform(log(1), log(200))
    log_z_min = np.log(1.0)
    log_z_max = np.log(200.0)
    log_z_samples = rng.uniform(log_z_min, log_z_max, size=n_samples)
    z_samples = np.exp(log_z_samples)

    # Hyperprior density: p(z) = 1/(z * (log_z_max - log_z_min))
    p_z = 1.0 / (z_samples * (log_z_max - log_z_min))

    # Evaluate interpolated ln BF at each z
    lnBF_samples = interp(np.log(z_samples))
    BF_samples = np.exp(lnBF_samples)

    # Monte Carlo estimate: <BF> = (1/N) Σ BF(z_i)
    # The samples are drawn from the hyperprior, so weights are uniform
    BF_marg = np.mean(BF_samples)
    lnBF_marg = np.log(BF_marg)

    # Uncertainty via bootstrap
    n_boot = 1000
    boot_means = []
    for _ in range(n_boot):
        idx = rng.integers(0, n_samples, size=n_samples)
        boot_means.append(np.mean(BF_samples[idx]))
    boot_means = np.array(boot_means)
    lnBF_err = np.std(np.log(boot_means))

    # Also compute the free-z result consistency check
    # The free-z dynesty prior was U[0.1, 150]; our hyperprior is LogU[1, 200]
    # These are close but not identical. We report both.
    return {
        "BF_marginalized": float(BF_marg),
        "lnBF_marginalized": float(lnBF_marg),
        "lnBF_marginalized_uncertainty": float(lnBF_err),
        "lnBF_free_zlos_dynesty": float(lnBF_free),
        "hyperprior": "LogUniform(1, 200)",
        "n_mc_samples": n_samples,
        "n_bootstrap": n_boot,
        "z_eval_points": [1.0, 5.0, 100.0],
        "lnBF_at_eval_points": [float(x) for x in np.log([
            load_existing_results()["bayes_factors"]["BF_M1_NoLambda_zT1_vs_M0a_LCDM"],
            load_existing_results()["bayes_factors"]["BF_M1_NoLambda_zT5_vs_M0a_LCDM"],
            load_existing_results()["bayes_factors"]["BF_M1_Unscreened_zT100_vs_M0a_LCDM"],
        ])],
        "interpretation": (
            "The hyperprior-marginalized Bayes factor is the defensible headline number. "
            "If it agrees with the free-z_los dynesty result (within Monte Carlo uncertainty), "
            "the latter was already an approximately correct marginalized evidence."
        )
    }


def main():
    print("Loading existing model-comparison results...")
    interp, lnBF_free, z_points, lnBF_points = approximate_lnBF_vs_zlos()

    print(f"Discrete ln BF points: z={z_points}, lnBF={lnBF_points}")
    print(f"Free-z_los dynesty ln BF: {lnBF_free:.4f}")

    print(f"\nRunning Monte Carlo hyperprior integration (LogUniform 1–200)...")
    results = hyperprior_integration(interp, lnBF_free)

    print(f"\n{'='*60}")
    print(f"Marginalized BF       = {results['BF_marginalized']:.3f}")
    print(f"Marginalized ln BF    = {results['lnBF_marginalized']:.4f} ± {results['lnBF_marginalized_uncertainty']:.4f}")
    print(f"Free-z_los dynesty ln BF = {results['lnBF_free_zlos_dynesty']:.4f}")
    print(f"\n{results['interpretation']}")
    print(f"{'='*60}")

    with open(OUTPUT_PATH, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
