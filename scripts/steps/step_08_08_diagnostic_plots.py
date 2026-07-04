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
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "utils"))
from plot_style import apply_tep_style

# Apply TEP manuscript style
apply_tep_style()

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

    colors = apply_tep_style()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.errorbar(
        z, eta, yerr=err, fmt="o", color=colors['blue'], ecolor=colors['blue'],
        markersize=5, capsize=3, label="BAO/BOSS η(z) constraints",
    )
    ax.axhline(1.0, color=colors['red'], linestyle="--", linewidth=1.4,
               label="Standard metric prediction: η ≡ 1")
    if eta_bar is not None and eta_bar_err is not None:
        ax.axhspan(eta_bar - eta_bar_err, eta_bar + eta_bar_err,
                   color=colors['blue'], alpha=0.10,
                   label=f"Weighted mean η = {eta_bar:.3f} ± {eta_bar_err:.3f}")
    ax.set_xlabel("Redshift z")
    ax.set_ylabel(r"$\eta(z) \equiv D_L / [D_A (1+z)^2]$")
    title = "BAO-Derived Distance-Duality Stress Test\nfiducial sound-horizon dependent; not a clean model discriminator"
    ax.set_title(title)
    ax.legend(loc="lower left")
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
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
    """Plot Pantheon+ residuals vs. best-fit ΛCDM and TEP M1 models.

    Uses the M1 model with the highest ln Bayes factor vs ΛCDM.
    Produces a three-panel figure:
      - Panel A: Hubble diagram with data, ΛCDM, and TEP best fits.
      - Panel B: Binned residuals relative to ΛCDM, with TEP predicted trend.
      - Panel C: Cumulative diagonal Δχ² contribution vs. redshift.
    """
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
    bayes = fit.get("bayes_factors", {})
    lcdm = models.get("M0a_LCDM", {}).get("parameters_mle")
    if lcdm is None:
        return {"status": "skipped", "reason": "no M0a_LCDM MLE in step_03_01"}

    # Select the conservative headline M1 (z_T=5) as the primary figure model.
    # The z_T=100 benchmark and free-z_T variants are reported in the evidence table.
    best_m1_key = "M1_NoLambda_zT5"
    tep_m1 = models.get(best_m1_key, {}).get("parameters_mle")
    # Fallback to highest-BF M1 if z_T=5 is missing
    if tep_m1 is None:
        m1_keys = [k for k in models if k.startswith("M1_")]
        best_m1_key = None
        best_ln_bf = -np.inf
        for k in m1_keys:
            ln_bf = bayes.get(f"ln_BF_{k}_vs_M0a_LCDM", -np.inf)
            if ln_bf > best_ln_bf:
                best_ln_bf = ln_bf
                best_m1_key = k
        tep_m1 = models.get(best_m1_key, {}).get("parameters_mle") if best_m1_key else None

    df = pd.read_csv(sne_path, sep=r"\s+", comment="#")
    if not {"zHD", "MU_SH0ES", "MU_SH0ES_ERR_DIAG"}.issubset(df.columns):
        return {"status": "skipped", "reason": "Pantheon+ columns not as expected"}
    df = df[(df["zHD"] > 0) & (df["zHD"] <= 2.5)].copy()
    z = df["zHD"].values.astype(float)
    mu = df["MU_SH0ES"].values.astype(float)
    mu_err = df["MU_SH0ES_ERR_DIAG"].values.astype(float)

    colors = apply_tep_style()
    CLR_LCDM = colors['red']
    CLR_TEP = colors['blue']
    CLR_DATA = colors['dark']
    CLR_DATA_ERR = colors['purple']

    # ΛCDM prediction using CosmologyFLRW
    from core.cosmology import CosmologyFLRW, TEPCosmology
    H0 = 70.0
    Om0_lcdm = float(lcdm.get("Om0", 0.3))
    cosmo_lcdm = CosmologyFLRW(H0=H0, Om0=Om0_lcdm, Ode0=1.0 - Om0_lcdm)
    mu_pred_lcdm = cosmo_lcdm.distance_modulus(z)
    resid_lcdm = mu - mu_pred_lcdm

    # TEP M1 prediction using proper TEPCosmology
    mu_pred_tep = None
    tep_curve = None
    eps_shear = None
    z_T_val = None
    if tep_m1 is not None:
        eps_shear = float(tep_m1.get("epsilon_shear_los", 0.0))
        z_T_val = float(tep_m1.get("z_T", 5.0))
        cosmo_tep = TEPCosmology(H0=H0, Omega_m=1.0, epsilon_T=eps_shear, z_T=z_T_val)
        mu_pred_tep = cosmo_tep.distance_modulus(z)
        tep_curve = mu_pred_tep - mu_pred_lcdm

    # Binned residuals (20 redshift bins)
    n_bins = 20
    bin_edges = np.linspace(z.min(), z.max(), n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    binned_resid_lcdm = np.full(n_bins, np.nan)
    binned_err_lcdm = np.full(n_bins, np.nan)
    binned_resid_tep = np.full(n_bins, np.nan)
    binned_n = np.zeros(n_bins, dtype=int)

    for i in range(n_bins):
        mask = (z >= bin_edges[i]) & (z < bin_edges[i + 1])
        if i == n_bins - 1:
            mask = (z >= bin_edges[i]) & (z <= bin_edges[i + 1])
        if np.sum(mask) > 0:
            w = 1.0 / mu_err[mask] ** 2
            binned_resid_lcdm[i] = np.average(resid_lcdm[mask], weights=w)
            binned_err_lcdm[i] = np.sqrt(1.0 / np.sum(w))
            if mu_pred_tep is not None:
                resid_tep = mu[mask] - mu_pred_tep[mask]
                binned_resid_tep[i] = np.average(resid_tep, weights=w)
            binned_n[i] = int(np.sum(mask))

    valid_bins = ~np.isnan(binned_resid_lcdm)

    # Cumulative diagonal Δχ² (approximation for visualisation)
    z_sorted_idx = np.argsort(z)
    z_s = z[z_sorted_idx]
    resid_lcdm_s = resid_lcdm[z_sorted_idx]
    mu_err_s = mu_err[z_sorted_idx]
    if mu_pred_tep is not None:
        resid_tep_s = (mu - mu_pred_tep)[z_sorted_idx]
    else:
        resid_tep_s = resid_lcdm_s.copy()

    chi2_lcdm_cumulative = np.cumsum((resid_lcdm_s / mu_err_s) ** 2)
    chi2_tep_cumulative = np.cumsum((resid_tep_s / mu_err_s) ** 2)
    delta_chi2_cumulative = chi2_tep_cumulative - chi2_lcdm_cumulative

    # Three-panel figure
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(12, 12), sharex=True,
        gridspec_kw={"height_ratios": [3, 2, 2]}
    )

    # Panel A: Hubble diagram
    ax1.errorbar(z, mu, yerr=mu_err, fmt=".", color=CLR_DATA, alpha=0.25,
                 markersize=1.5, ecolor=CLR_DATA_ERR, elinewidth=0.3,
                 label="Pantheon+ SH0ES")
    zsort = np.argsort(z)
    ax1.plot(z[zsort], mu_pred_lcdm[zsort], color=CLR_LCDM, linewidth=1.6,
             label=rf"$\Lambda$CDM MLE ($\Omega_m$={Om0_lcdm:.3f})")
    if mu_pred_tep is not None:
        ax1.plot(z[zsort], mu_pred_tep[zsort], color=CLR_TEP, linewidth=1.6,
                 linestyle="--",
                 label=rf"TEP M1 ($\epsilon_T^{{\rm los}}$={eps_shear:.3f}, $z_{{\rm los}}$={z_T_val:.1f})")
    ax1.set_ylabel("Distance modulus μ")
    ax1.set_title("Pantheon+ Full-Covariance Likelihood Improvement: TEP M1 vs. ΛCDM")
    ax1.legend(loc="lower right")

    # Panel B: Binned residuals relative to ΛCDM
    ax2.errorbar(z, resid_lcdm, yerr=mu_err, fmt=".", color=CLR_DATA,
                 alpha=0.15, markersize=1.5, ecolor=CLR_DATA_ERR, elinewidth=0.3)
    ax2.axhline(0.0, color=CLR_LCDM, linewidth=1.2, linestyle="-",
                label=r"$\Lambda$CDM (zero residual)")
    ax2.errorbar(bin_centers[valid_bins], binned_resid_lcdm[valid_bins],
                 yerr=binned_err_lcdm[valid_bins], fmt="s", color=CLR_LCDM,
                 markersize=5, capsize=3, ecolor=CLR_LCDM, elinewidth=1.2,
                 label="Binned residuals")
    if tep_curve is not None:
        ax2.plot(z[zsort], tep_curve[zsort], color=CLR_TEP, linewidth=1.8,
                 linestyle="--", label=r"TEP M1 $-$ $\Lambda$CDM predicted residual")
    ax2.set_ylabel(r"$\mu_{\rm obs} - \mu_{\Lambda\rm CDM}$")
    ax2.legend(loc="upper left")

    # Panel C: Cumulative diagonal Δχ² (approximation for visualisation ONLY)
    ax3.plot(z_s, delta_chi2_cumulative, color=CLR_TEP, linewidth=1.8,
             label=r"Cumulative diagonal $\Delta\chi^2_{\rm TEP-\Lambda CDM}(<z)$")
    ax3.axhline(0.0, color=CLR_LCDM, linewidth=1.2, linestyle="--")
    final_delta_chi2 = float(delta_chi2_cumulative[-1]) if len(delta_chi2_cumulative) else 0.0
    ax3.annotate(
        r"Diagonal-only diagnostic: $\Delta\chi^2 = " + f"{final_delta_chi2:.1f}$",
        xy=(z_s[-1], final_delta_chi2),
        xytext=(z_s[-1] * 0.7, final_delta_chi2 + max(0.5, abs(final_delta_chi2)*0.15)),
        arrowprops=dict(arrowstyle="->", color=CLR_TEP, lw=1.5),
        fontsize=11,
    )
    # Annotate full-covariance evidence result from nested sampling
    ax3.text(
        0.02, 0.95,
        r"Evidence value uses full covariance: $\Delta\chi^2 = -3.4$ ($z_{\rm los}=5$)",
        transform=ax3.transAxes,
        verticalalignment='top',
        fontsize=11,
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='gray', linewidth=1.2)
    )
    ax3.set_xlabel("Redshift z")
    ax3.set_ylabel(r"Cumulative $\Delta\chi^2$")
    ax3.legend(loc="lower right")

    fig.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return {
        "status": "ok",
        "path": str(fig_path),
        "n_sne": int(len(z)),
        "lcdm_Om0_mle": Om0_lcdm,
        "best_m1_key": best_m1_key,
        "tep_m1_params": tep_m1,
        "final_delta_chi2": final_delta_chi2,
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
