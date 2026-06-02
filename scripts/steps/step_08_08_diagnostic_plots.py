#!/usr/bin/env python3
"""Step 088: Diagnostic Plots generation for Manuscript.

Produces real, data-driven diagnostic figures from upstream pipeline outputs:
- ``distance_duality.png`` : observational η(z) constraints from
  ``step_04_04_distance_duality.json`` overlaid with the TEP and ΛCDM
  predictions (both analytically η=1 by metric construction).
- ``hubble_residuals.png`` : Pantheon+ μ residuals against the best-fit
  ΛCDM model from ``step_03_01_three_model_comparison.json`` with the TEP
  M1 (no Lambda) residuals overlaid. Skipped if Pantheon+ data missing.

Any figure whose required upstream artefact is missing is *skipped*
(reason recorded in the payload) rather than replaced with a synthetic image.
Fabrication is strictly prohibited.
"""

import json
import os
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from c0_common import (
    RESULTS_DIR,
    TEPLogger,
    ensure_dirs,
    print_status,
    set_step_logger,
    step_json_path,
    write_json,
)

STEP_ID = "step_08_08_diagnostic_plots"
FIG_DIR = Path("results/figures")


def _plot_distance_duality(fig_path: Path) -> dict:
    """Plot observational distance-duality constraints from step_04_04."""
    src = RESULTS_DIR / "step_04_04_distance_duality.json"
    if not src.exists():
        return {"status": "skipped", "reason": f"{src} not found"}

    with open(src) as f:
        ddr = json.load(f)

    constraints = ddr.get("ddr_constraints", [])
    if not constraints:
        return {"status": "skipped", "reason": "no ddr_constraints in step_04_04"}

    z = np.array([c["z"] for c in constraints])
    eta = np.array([c["eta_obs"] for c in constraints])
    err = np.array([c["eta_err"] for c in constraints])

    stats = ddr.get("statistical_results", {})
    eta_bar = stats.get("eta_weighted_mean")
    eta_bar_err = stats.get("eta_weighted_err")
    dev_sigma = stats.get("deviation_from_unity_sigma")

    fig, ax = plt.subplots(figsize=(8.5, 5.0))
    ax.errorbar(
        z, eta, yerr=err, fmt="o", color="#1f4e8c", ecolor="#1f4e8c",
        markersize=5, capsize=3, label="BAO/BOSS η(z) constraints",
    )
    ax.axhline(1.0, color="#b21f1f", linestyle="--", linewidth=1.4,
               label="TEP / ΛCDM prediction: η ≡ 1 (metric-preserved)")
    if eta_bar is not None and eta_bar_err is not None:
        ax.axhspan(eta_bar - eta_bar_err, eta_bar + eta_bar_err,
                   color="#1f4e8c", alpha=0.10,
                   label=f"Weighted mean η = {eta_bar:.3f} ± {eta_bar_err:.3f}")
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$\eta(z) \equiv D_L / [D_A (1+z)^2]$")
    title = "Distance-Duality Relation: Observations vs. TEP/ΛCDM"
    if dev_sigma is not None:
        title += f"  (data: {dev_sigma:.1f}σ from unity)"
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower left", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "ok",
        "path": str(fig_path),
        "n_points": int(len(z)),
        "eta_weighted_mean": eta_bar,
        "eta_weighted_err": eta_bar_err,
        "deviation_from_unity_sigma": dev_sigma,
        "source": str(src),
    }


def _plot_hubble_residuals(fig_path: Path) -> dict:
    """Plot Pantheon+ residuals vs. best-fit ΛCDM and TEP M1 models."""
    sne_path = Path("data/raw/pantheon_plus_shoes.dat")
    fit_path = RESULTS_DIR / "step_03_01_three_model_comparison.json"
    if not sne_path.exists():
        return {"status": "skipped", "reason": f"{sne_path} not found"}
    if not fit_path.exists():
        return {"status": "skipped", "reason": f"{fit_path} not found"}

    try:
        import pandas as pd
    except ImportError:
        return {"status": "skipped", "reason": "pandas not available"}

    with open(fit_path) as f:
        fit = json.load(f)
    models = fit.get("models", {})
    lcdm = models.get("M0a_LCDM", {}).get("parameters_mle")
    tep_m1 = (
        models.get("M1_free_zT", {}).get("parameters_mle")
        or models.get("M1_NoLambda_zT5", {}).get("parameters_mle")
        or models.get("M1_NoLambda_zT1", {}).get("parameters_mle")
    )
    if lcdm is None:
        return {"status": "skipped", "reason": "no M0a_LCDM MLE in step_03_01"}

    df = pd.read_csv(sne_path, sep=r"\s+", comment="#")
    if not {"zHD", "MU_SH0ES", "MU_SH0ES_ERR_DIAG"}.issubset(df.columns):
        return {"status": "skipped", "reason": "Pantheon+ columns not as expected"}
    df = df[(df["zHD"] > 0) & (df["zHD"] <= 2.5)].copy()
    z = df["zHD"].values
    mu = df["MU_SH0ES"].values
    mu_err = df["MU_SH0ES_ERR_DIAG"].values

    # Predicted apparent magnitude m_b from ΛCDM MLE: μ = m_b - M
    # The MLE stores Om0 and M; we compute μ_pred via flat ΛCDM luminosity distance.
    from scipy.integrate import cumulative_trapezoid
    c_km_s = 2.99792458e5
    H0 = 70.0  # MLE fit is performed at fixed H0 = 70 (see step_03_01)
    Om0_lcdm = float(lcdm.get("Om0", 0.3))

    def mu_flat_lcdm(zv, Om0):
        zg = np.linspace(0, float(np.max(zv)) + 1e-3, 4096)
        Ez = np.sqrt(Om0 * (1 + zg) ** 3 + (1 - Om0))
        DC = cumulative_trapezoid(c_km_s / (H0 * Ez), zg, initial=0.0)
        DL = (1 + zv) * np.interp(zv, zg, DC)
        return 5.0 * np.log10(np.maximum(DL, 1e-12)) + 25.0

    mu_pred_lcdm = mu_flat_lcdm(z, Om0_lcdm)
    resid_lcdm = mu - mu_pred_lcdm

    tep_curve = None
    if tep_m1 is not None:
        Om0_tep = float(tep_m1.get("Om0", 1.0))  # M1_NoLambda has Om0 fixed to 1.0
        eps_T = float(tep_m1.get("epsilon_T", 0.0))
        z_T = float(tep_m1.get("z_T", 1.0))
        # M1 prediction: ΛCDM-flat with Om0_tep plus TEP transport term
        # μ_TEP(z) = μ_flat(z; Om0_tep) + (5/ln10) * ε_T * (1 - exp(-z/z_T))
        mu_pred_tep = mu_flat_lcdm(z, Om0_tep) + (5.0 / np.log(10.0)) * eps_T * (1.0 - np.exp(-z / max(z_T, 1e-6)))
        tep_curve = mu_pred_tep - mu_pred_lcdm  # residual relative to LCDM

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.0, 7.0), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2]})

    ax1.errorbar(z, mu, yerr=mu_err, fmt=".", color="#444", alpha=0.35,
                 markersize=2, ecolor="#888", elinewidth=0.4, label="Pantheon+ SH0ES")
    zsort = np.argsort(z)
    ax1.plot(z[zsort], mu_pred_lcdm[zsort], color="#b21f1f", linewidth=1.6,
             label=f"ΛCDM MLE (Ωₘ={Om0_lcdm:.3f})")
    if tep_curve is not None:
        ax1.plot(z[zsort], (mu_pred_lcdm + tep_curve)[zsort], color="#1f4e8c",
                 linewidth=1.6, linestyle="--",
                 label=f"TEP M1 (εₜ={eps_T:.3f}, z_T={z_T:.2f})")
    ax1.set_ylabel("Distance modulus μ")
    ax1.set_title("Pantheon+ Hubble Diagram and Residuals vs. ΛCDM Best Fit")
    ax1.legend(loc="lower right", framealpha=0.95)
    ax1.grid(True, alpha=0.3)

    ax2.errorbar(z, resid_lcdm, yerr=mu_err, fmt=".", color="#444",
                 alpha=0.35, markersize=2, ecolor="#888", elinewidth=0.4)
    ax2.axhline(0.0, color="#b21f1f", linewidth=1.2)
    if tep_curve is not None:
        ax2.plot(z[zsort], tep_curve[zsort], color="#1f4e8c", linewidth=1.6,
                 linestyle="--", label="TEP M1 − ΛCDM")
        ax2.legend(loc="upper left", framealpha=0.95)
    ax2.set_xlabel("Redshift z")
    ax2.set_ylabel("μ − μ_ΛCDM")
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "ok",
        "path": str(fig_path),
        "n_sne": int(len(z)),
        "lcdm_Om0_mle": Om0_lcdm,
        "tep_m1_params": tep_m1,
    }


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    figures = {}

    print_status("Generating distance-duality figure...", "PROCESS")
    try:
        figures["distance_duality"] = _plot_distance_duality(FIG_DIR / "distance_duality.png")
    except Exception as exc:
        figures["distance_duality"] = {"status": "error", "reason": str(exc)}

    print_status("Generating Hubble residuals figure...", "PROCESS")
    try:
        figures["hubble_residuals"] = _plot_hubble_residuals(FIG_DIR / "hubble_residuals.png")
    except Exception as exc:
        figures["hubble_residuals"] = {"status": "error", "reason": str(exc)}

    generated = [v["path"] for v in figures.values() if v.get("status") == "ok"]
    skipped = {k: v for k, v in figures.items() if v.get("status") != "ok"}

    payload = {
        "step": STEP_ID,
        "status": "completed" if generated else "blocked",
        "generated_plots": generated,
        "figures": figures,
        "skipped": skipped,
        "validation": {
            "claim_gate": "open" if not skipped else "partial",
            "research_grade": bool(generated and not skipped),
            "note": (
                "Diagnostic plots are generated only from real upstream artefacts. "
                "Synthetic image figures are forbidden."
            ),
        },
        "timestamp": int(time.time()),
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(
        f"Step {STEP_ID}: {len(generated)} figure(s) generated, {len(skipped)} skipped",
        "SUCCESS" if generated else "WARNING",
    )
    return payload


if __name__ == "__main__":
    run()
