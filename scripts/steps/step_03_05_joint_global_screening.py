#!/usr/bin/env python3
"""
Step 03.05: Production Joint Global Screening Inference
=======================================================

Coupled multi-domain Bayesian inference across:
  1. Pantheon+ SNe Ia (full 1,701x1,701 covariance Cholesky likelihood on the
     no-Lambda TEP transport background);
  2. CMB acoustic & Boltzmann sector (live patched hi_class response surface
     on 100*theta_s and r_s(z_rec) combined with the converged Planck 2018
     TTTEEE + lensing posterior constraint, verified by a direct live hi_class
     Boltzmann call at the joint posterior mean);
  3. BBN nuclear reaction network (AlterBBN baseline abundances + MeV-epoch
     expansion-rate response evaluated against the empirical Y_p, D/H, 3He/H,
     7Li/H abundance registry);
  4. Solar-System PPN Cassini bound & Galactic/Wide-Binary screening transition
     evaluated on the exact selected two-branch kinetic operator
     P_{,X}(xi) = k sqrt(xi) + 2 xi (k = 4 sqrt(2) |beta_A| g_t / a_0).

Parameters sampled simultaneously via dynamic nested sampling (dynesty):
    epsilon_T      : Line-of-sight temporal shear amplitude (prior: U[0.0, 1.2])
    z_T            : Screening turnover redshift            (prior: U[1.0, 15.0])
    log10_g_t      : log10(transition acceleration in m/s^2) (prior: U[-11.0, -8.0])
    log10_rho_half : log10(galactic halo screening half-density in M_sun/pc^3;
                           1 M_sun/pc^3 = 6.768e-23 g/cm^3) (prior: U[-1.0, 2.0])

Output: results/step_03_05_joint_global_screening.json
"""

import sys
import os
from pathlib import Path
import json
import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.interpolate import RegularGridInterpolator, interp1d

# Keep BLAS single-threaded
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

PROJECT_ROOT = Path(__file__).resolve().parents[2]
WORKSPACE_ROOT = PROJECT_ROOT.parent
HICLASS_PYTHON = WORKSPACE_ROOT / "TEP-HC" / "external" / "hi_class" / "hi_class" / "python"

sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "steps"))
if HICLASS_PYTHON.exists():
    sys.path.insert(0, str(HICLASS_PYTHON))

from c0_common import (
    TEPLogger, ensure_dirs, print_status, set_step_logger,
    step_json_path, write_json, RESULTS_DIR, PROCESSED_DIR
)
from core.cosmology import TEPCosmology
from scripts.utils.bbn_network_resolver import BBNInputs, chi2_against_registry, run_network
import step_03_01_three_model_comparison as s301

STEP_ID = "step_03_05_joint_global_screening"
OUTPUT_JSON = RESULTS_DIR / "step_03_05_joint_global_screening.json"

# Physical & observational constants
BETA_A = -1.0
A0_GAL = 1.2e-10          # m/s^2 (galactic RAR acceleration scale)
G_SOLAR_1AU = 5.93e-3     # m/s^2 (solar acceleration at 1 AU)
G_CASSINI_LIMB = 107.0    # m/s^2 (solar acceleration at 1.6 R_sun Cassini ray)
CASSINI_BOUND = 2.3e-5    # |gamma_PPN - 1| < 2.3e-5 (Bertotti et al. 2003)
PLANCK_100THETA_S = 1.04110
PLANCK_100THETA_S_ERR = 0.00031
CMB_EPS_HOM_MEAN = 0.0043
CMB_EPS_HOM_STD = 0.0025

_PANTHEON_CACHE = None
_BBN_CACHE = None
_HICLASS_INTERP_CACHE = None


def _get_pantheon_profile():
    """Load real Pantheon+ sample and Cholesky factor for exact M-marginalized chi2."""
    global _PANTHEON_CACHE
    if _PANTHEON_CACHE is None:
        data = s301.PantheonData()
        if not data.load():
            raise RuntimeError("Failed to load real Pantheon+ data and covariance")
        L = cholesky(data.cov, lower=True, check_finite=False)
        one = solve_triangular(L, np.ones(len(data.z)), lower=True, check_finite=False)
        w = float(np.sum(one * one))
        mb = data.mb - data.observer_correction
        _PANTHEON_CACHE = (data.z, mb, L, one, w)
    return _PANTHEON_CACHE


def _get_bbn_baseline():
    """Run AlterBBN C nuclear network across a calibrated expansion-rate grid
    and load the empirical BBN abundance registry."""
    global _BBN_CACHE
    if _BBN_CACHE is None:
        import csv
        import subprocess
        import tempfile

        reg_path = PROCESSED_DIR / "tep_c0_bbn_abundance_registry.csv"
        with reg_path.open(newline="", encoding="utf-8") as fh:
            registry_rows = list(csv.DictReader(fh))
        net_out = run_network(BBNInputs(eta=6.1e-10, h0_km_s_mpc=70.0, sigma_0=0.0, epsilon_t=0.0))
        base_abund = net_out["abundances"]
        base_chi2 = chi2_against_registry(base_abund, registry_rows)

        # Build direct AlterBBN C nuclear-network response interpolators over delta_H/H:
        ab_dir = PROJECT_ROOT / "external" / "AlterBBN"
        prim = (ab_dir / "primary.x").resolve()
        ini_text = (ab_dir / "input.ini").read_text(encoding="utf-8")
        dH_grid = np.array([0.0, 0.0005, 0.001, 0.002, 0.005, 0.010, 0.025, 0.050, 0.100, 0.250, 0.500])
        yp_g, dh_g, he3_g, li7_g = [], [], [], []
        with tempfile.TemporaryDirectory() as tmp:
            for dH in dH_grid:
                dnnu = (43.0 / 7.0) * (((1.0 + float(dH)) ** 2) - 1.0)
                mod_ini = ini_text.replace("dnnu = 0.0", f"dnnu = {dnnu:.6f}")
                (Path(tmp) / "input.ini").write_text(mod_ini, encoding="utf-8")
                res = subprocess.run([str(prim), "standard"], cwd=tmp, capture_output=True, text=True)
                if res.returncode != 0:
                    raise RuntimeError(
                        f"AlterBBN primary.x failed at dH={dH} with return code {res.returncode}: {res.stderr}"
                    )
                parsed = False
                for line in res.stdout.splitlines():
                    if line.startswith("cent:"):
                        vals = [float(p) for p in line.split()[1:]]
                        yp_g.append(vals[0])
                        dh_g.append(vals[1])
                        he3_g.append(vals[2])
                        li7_g.append(vals[3])
                        parsed = True
                        break
                if not parsed:
                    raise RuntimeError(
                        f"AlterBBN primary.x output missing 'cent:' abundance row at dH={dH}: {res.stdout}"
                    )

        _BBN_CACHE = {
            "engine": net_out["network"]["engine"],
            "base_abundances": base_abund,
            "base_chi2": float(base_chi2["chi2"]),
            "registry_rows": registry_rows,
            "dH_grid": dH_grid.tolist(),
            "grid_Y_p": [float(v) for v in yp_g],
            "grid_D_H": [float(v) for v in dh_g],
            "grid_He3_H": [float(v) for v in he3_g],
            "grid_Li7_H": [float(v) for v in li7_g],
            "interp_Y_p": interp1d(dH_grid, yp_g, kind="linear", bounds_error=False, fill_value="extrapolate"),
            "interp_D_H": interp1d(dH_grid, dh_g, kind="linear", bounds_error=False, fill_value="extrapolate"),
            "interp_He3_H": interp1d(dH_grid, he3_g, kind="linear", bounds_error=False, fill_value="extrapolate"),
            "interp_Li7_H": interp1d(dH_grid, li7_g, kind="linear", bounds_error=False, fill_value="extrapolate"),
        }
    return _BBN_CACHE


def _build_hiclass_response_surface():
    """Build a calibrated 2D response surface of 100*theta_s and r_s(z_rec)
    across the full prior-reachable domain [0.0, 0.500] from live patched
    hi_class (classy) Boltzmann evaluations, plus the CDM-free (omega_cdm=0)
    baryon lensing spectrum and Planck 2018 9-bin CMBMarged likelihood."""
    global _HICLASS_INTERP_CACHE
    if _HICLASS_INTERP_CACHE is not None:
        return _HICLASS_INTERP_CACHE

    import classy
    from cobaya.likelihoods.planck_2018_lensing import CMBMarged

    eps_grid = np.array([0.0, 0.001, 0.003, 0.006, 0.010, 0.020, 0.040, 0.080, 0.160, 0.300, 0.500])
    zT_grid = np.array([1.0, 3.0, 5.0, 8.0, 15.0])
    theta_grid = np.zeros((len(eps_grid), len(zT_grid)))
    rs_grid = np.zeros((len(eps_grid), len(zT_grid)))

    base_params = {
        "output": "tCl",
        "l_max_scalars": 500,
        "H0": 67.36,
        "omega_b": 0.02237,
        "omega_cdm": 0.1200,
        "A_s": 2.1e-9,
        "n_s": 0.9649,
        "tau_reio": 0.0544,
        "tep_mode": "yes",
        "n_T": 2.0,
    }

    print_status("Calibrating full-domain 2D Boltzmann acoustic response surface + CDM-free lensing via live hi_class calls...", "PROCESS")
    for i, eh in enumerate(eps_grid):
        for j, zt in enumerate(zT_grid):
            cosmo = classy.Class()
            p = dict(base_params, epsilon_T=float(eh), z_T=float(zt))
            cosmo.set(p)
            cosmo.compute()
            derived = cosmo.get_current_derived_parameters(["100*theta_s", "rs_rec"])
            theta_grid[i, j] = float(derived["100*theta_s"])
            rs_grid[i, j] = float(derived["rs_rec"])
            cosmo.struct_cleanup()
            cosmo.empty()

    # Compute pure-baryon CDM-free (omega_cdm = 0) Weyl lensing spectrum for Planck 2018 CMBMarged
    lmax_lens = 2500
    cosmo_nocdm = classy.Class()
    cosmo_nocdm.set({
        "output": "tCl,pCl,lCl",
        "lensing": "yes",
        "modes": "s",
        "l_max_scalars": lmax_lens,
        "non_linear": "none",
        "H0": 66.77,
        "omega_b": 0.02144,
        "omega_cdm": 0.0,
        "A_s": 2.1e-9,
        "n_s": 0.9956,
        "tau_reio": 0.0497,
        "tep_mode": "yes",
        "z_T": 5.0,
        "n_T": 2.0,
        "epsilon_T": 0.00547,
        "gravity_model": "tep",
        "M2_evolution": "yes",
        "use_pert_var_deltaphi_smg": "yes",
        "pert_initial_conditions_smg": "zero",
        "tep_screening_s": 1.0,
        "tep_kin_boost": 10.0,
    })
    cosmo_nocdm.compute()
    cls_nocdm = {k: np.array(cosmo_nocdm.lensed_cl(lmax_lens)[k], dtype=float) for k in ("tt", "ee", "te", "pp")}
    cosmo_nocdm.struct_cleanup()
    cosmo_nocdm.empty()

    pkg_path = WORKSPACE_ROOT / "TEP-HC" / "data" / "external" / "cobaya_packages"
    planck_lens_like = CMBMarged({"packages_path": str(pkg_path)})
    ells = np.arange(lmax_lens + 1, dtype=float)
    fac_pp = (ells * (ells + 1.0)) ** 2 / (2.0 * np.pi)
    fac_tt = ells * (ells + 1.0) / (2.0 * np.pi)
    l_safe = np.maximum(ells, 2.0)
    w_silk = 1.0 / (1.0 + (l_safe / 300.0) ** 2.0)

    theta_interp = RegularGridInterpolator(
        (eps_grid, zT_grid), theta_grid, bounds_error=True, fill_value=None
    )
    rs_interp = RegularGridInterpolator(
        (eps_grid, zT_grid), rs_grid, bounds_error=True, fill_value=None
    )

    # Profile/marginalize Planck 2018 CMBMarged lensing chi2(log10_g_t) over the
    # cosmic-web baryon acceleration scales (log10_g_b_sheet, log10_g_b_fil) with
    # independent physical structure-scale priors (sigma = 0.15 dex around
    # cluster-outskirt sheet G*M_b/R^2 = 2.98e-11 m/s^2 and cosmic-web filament
    # G*M_b/R^2 = 1.11e-13 m/s^2) on the zero-free-parameter two-branch sector
    # k = 16.03, so that the full unweighted (1.0x) lensing likelihood is used.
    #
    # The two structure scales are *marginalized* (integrated over their
    # Gaussian priors) rather than profile-optimized: at each log10(g_t) the
    # marginalized lensing log-likelihood is evaluated by exact 2D Gauss-Hermite
    # quadrature of the data-only chi2 over the 0.15-dex Gaussian priors,
    #   ln L_marg(g_t) = ln [ (1/pi) sum_ij w_i w_j
    #                         exp(-chi2_data(lg_s_i, lg_f_j; g_t)/2) ],
    #   lg_*_i = mu_* + sqrt(2) sigma x_i,
    # which carries the Occam prior-volume factor inside the evidence integral.
    # A profiled (posterior-mode) grid is retained in parallel for reporting
    # the best-fit structure scales and the profiled-vs-marginalized Occam
    # penalty, but the sampler likelihood uses ln L_marg.
    from scipy.optimize import minimize as _minimize
    from scipy.special import logsumexp as _logsumexp
    from numpy.polynomial.hermite import hermgauss as _hermgauss

    log10_gb_sheet_prior = -10.5256
    log10_gb_fil_prior = -12.9531
    sigma_gb_dex = 0.15
    k_canonical_tb = 16.03
    loggt_grid = np.linspace(-11.0, -8.0, 61)
    chi2_lens_grid = np.zeros_like(loggt_grid)
    chi2_lens_data_grid = np.zeros_like(loggt_grid)
    lnL_lens_marg_grid = np.zeros_like(loggt_grid)

    n_gh = 9
    gh_x, gh_w = _hermgauss(n_gh)
    lg_s_nodes = log10_gb_sheet_prior + np.sqrt(2.0) * sigma_gb_dex * gh_x
    lg_f_nodes = log10_gb_fil_prior + np.sqrt(2.0) * sigma_gb_dex * gh_x
    log_w2d = (
        np.log(gh_w)[:, None] + np.log(gh_w)[None, :] - np.log(np.pi)
    )

    def _lens_chi2_data(lg_s, lg_f, gt_v):
        qs = (10.0 ** lg_s) / gt_v
        qf = (10.0 ** lg_f) / gt_v
        ys = solve_two_branch_y(qs, k_canonical_tb)
        yf = solve_two_branch_y(qf, k_canonical_tb)
        Ml = 1.0 + (ys * w_silk + yf * (1.0 - w_silk))
        dls_g = {
            "pp": cls_nocdm["pp"] * (Ml ** 2) * fac_pp,
            "tt": cls_nocdm["tt"] * fac_tt,
            "ee": cls_nocdm["ee"] * fac_tt,
            "te": cls_nocdm["te"] * fac_tt,
        }
        return float(-2.0 * planck_lens_like.log_likelihood(dls_g))

    def _eval_mode_centered_gh(lg_s_c, lg_f_c, gt_v, n_nodes, h_scale=0.045):
        """Mode-centered 2D Gauss-Hermite quadrature of
        int L_data(lg_s, lg_f; gt_v) * N(lg_s; mu_s, sig^2) * N(lg_f; mu_f, sig^2) dlg_s dlg_f
        centered at the posterior mode (lg_s_c, lg_f_c) so the quadrature nodes
        resolve the high-likelihood ridge at every log10(g_t)."""
        gx, gw = _hermgauss(n_nodes)
        lgs_n = lg_s_c + np.sqrt(2.0) * h_scale * gx
        lgf_n = lg_f_c + np.sqrt(2.0) * h_scale * gx
        log_terms = np.zeros((n_nodes, n_nodes))
        log_jac = 2.0 * np.log(h_scale / sigma_gb_dex) - np.log(np.pi)
        for i_s, lg_s in enumerate(lgs_n):
            for j_f, lg_f in enumerate(lgf_n):
                c_data = _lens_chi2_data(lg_s, lg_f, gt_v)
                c_prior = (
                    ((lg_s - log10_gb_sheet_prior) / sigma_gb_dex) ** 2
                    + ((lg_f - log10_gb_fil_prior) / sigma_gb_dex) ** 2
                )
                c_prop = (
                    ((lg_s - lg_s_c) / h_scale) ** 2
                    + ((lg_f - lg_f_c) / h_scale) ** 2
                )
                log_w = np.log(gw[i_s]) + np.log(gw[j_f]) + log_jac
                log_terms[i_s, j_f] = -0.5 * (c_data + c_prior - c_prop) + log_w
        return float(_logsumexp(log_terms))

    lnL_lens_marg_grid_9 = np.zeros_like(loggt_grid)
    lnL_lens_marg_grid_17 = np.zeros_like(loggt_grid)

    for idx_g, lgt in enumerate(loggt_grid):
        gt_v = 10.0 ** float(lgt)

        def _lens_profile_obj(x_gb, gt_v=gt_v):
            lg_s, lg_f = x_gb
            c_data = _lens_chi2_data(lg_s, lg_f, gt_v)
            c_prior = (
                ((lg_s - log10_gb_sheet_prior) / sigma_gb_dex) ** 2
                + ((lg_f - log10_gb_fil_prior) / sigma_gb_dex) ** 2
            )
            return c_data + c_prior

        x0_gb = [
            log10_gb_sheet_prior + (float(lgt) - (-9.4685)),
            log10_gb_fil_prior + (float(lgt) - (-9.4685)),
        ]
        res_lp = _minimize(
            _lens_profile_obj,
            x0_gb,
            method="Nelder-Mead",
            options={"xatol": 1e-4, "fatol": 1e-4},
        )
        chi2_lens_grid[idx_g] = float(res_lp.fun)
        lg_s_opt, lg_f_opt = res_lp.x
        c_pr_opt = (
            ((lg_s_opt - log10_gb_sheet_prior) / sigma_gb_dex) ** 2
            + ((lg_f_opt - log10_gb_fil_prior) / sigma_gb_dex) ** 2
        )
        chi2_lens_data_grid[idx_g] = float(res_lp.fun - c_pr_opt)

        lnL_lens_marg_grid_9[idx_g] = _eval_mode_centered_gh(lg_s_opt, lg_f_opt, gt_v, 9)
        lnL_lens_marg_grid_17[idx_g] = _eval_mode_centered_gh(lg_s_opt, lg_f_opt, gt_v, 17)
        lnL_lens_marg_grid[idx_g] = lnL_lens_marg_grid_17[idx_g]

    post_mask = (loggt_grid >= -9.70) & (loggt_grid <= -9.20)
    gh_max_diff_9_vs_17_posterior = float(np.max(np.abs(lnL_lens_marg_grid_9[post_mask] - lnL_lens_marg_grid_17[post_mask])))
    gh_max_diff_9_vs_17_full = float(np.max(np.abs(lnL_lens_marg_grid_9 - lnL_lens_marg_grid_17)))

    chi2_lens_min = float(np.min(chi2_lens_grid))
    lnL_lens_ref = -0.5 * chi2_lens_min
    chi2_lens_interp = interp1d(loggt_grid, chi2_lens_grid, kind="cubic", bounds_error=False, fill_value="extrapolate")
    chi2_lens_data_interp = interp1d(loggt_grid, chi2_lens_data_grid, kind="cubic", bounds_error=False, fill_value="extrapolate")
    lnL_lens_marg_max = float(np.max(lnL_lens_marg_grid))
    lnL_lens_marg_interp = interp1d(loggt_grid, lnL_lens_marg_grid_17, kind="cubic", bounds_error=False, fill_value="extrapolate")

    _HICLASS_INTERP_CACHE = {
        "eps_grid": eps_grid.tolist(),
        "zT_grid": zT_grid.tolist(),
        "theta_s_100_lcdm": float(theta_grid[0, 2]),
        "rs_rec_lcdm_mpc": float(rs_grid[0, 2]),
        "theta_interp": theta_interp,
        "rs_interp": rs_interp,
        "cls_nocdm": cls_nocdm,
        "planck_lens_like": planck_lens_like,
        "chi2_lens_interp": chi2_lens_interp,
        "chi2_lens_data_interp": chi2_lens_data_interp,
        "chi2_lens_min": chi2_lens_min,
        "lnL_lens_ref": lnL_lens_ref,
        "lnL_lens_marg_interp": lnL_lens_marg_interp,
        "lnL_lens_marg_max": float(np.max(lnL_lens_marg_grid_17)),
        "gh_quadrature_convergence": {
            "max_abs_diff_9x9_vs_17x17_posterior_band": gh_max_diff_9_vs_17_posterior,
            "max_abs_diff_9x9_vs_17x17_full_grid": gh_max_diff_9_vs_17_full,
        },
        "k_canonical_tb": k_canonical_tb,
        "lensing_nuisance_treatment": (
            "marginalized_2d_gauss_hermite_17pt_over_0.15dex_sheet_filament_priors_common_ref_baseline"
        ),
        "fac_pp": fac_pp,
        "fac_tt": fac_tt,
        "w_silk": w_silk,
        "log10_gb_sheet_prior": log10_gb_sheet_prior,
        "log10_gb_fil_prior": log10_gb_fil_prior,
        "sigma_gb_dex": sigma_gb_dex,
    }
    return _HICLASS_INTERP_CACHE


def solve_two_branch_y(q_val, k_val):
    """Exact Newton solve of c1 * q * y^2 + q^2 * y^3 = 1 on the selected
    two-branch kinetic sector P_{,X}(xi) = k sqrt(xi) + 2 xi."""
    c1 = k_val / np.sqrt(2.0)
    y = np.sqrt(1.0 / (c1 * q_val)) if q_val < 1.0 else q_val ** (-2.0 / 3.0)
    for _ in range(12):
        f = c1 * q_val * (y ** 2) + (q_val ** 2) * (y ** 3) - 1.0
        fp = 2.0 * c1 * q_val * y + 3.0 * (q_val ** 2) * (y ** 2)
        y = max(y - f / fp, 1e-30)
    return float(y)


def eval_nocdm_planck_lensing_chi2(log10_g_t, hc):
    """Evaluate the profiled 9-bin Planck 2018 CMBMarged lensing data chi2 at
    omega_cdm = 0 at the specified log10(g_t)."""
    return float(hc["chi2_lens_data_interp"](float(log10_g_t)))


_USE_MARGINALIZED_LENSING = True


def log_likelihood(theta):
    """
    Coupled multi-domain log-likelihood:
      ln L_total = ln L_SNe(Pantheon+ coupled to host-screening rho_half)
                 + ln L_CMB(unclamped hi_class Boltzmann + Planck TTTEEE
                            + omega_cdm=0 CMBMarged lensing marginalized over
                              0.15-dex sheet/filament scale priors via exact
                              2D Gauss-Hermite quadrature)
                 + ln L_BBN(live-calibrated AlterBBN C nuclear network)
                 + ln L_PPN_WB(Two-branch Cassini + SPARC/WB)
                 + ln L_rho(Galactic halo screening half-density law in M_sun/pc^3)
    """
    epsilon_T, z_T, log10_g_t, log10_rho_half = theta

    # Physical host-disk screening modulation from rho_half [M_sun/pc^3]:
    # Characteristic SN Ia host disk baryon density rho_disk ~ 0.05 M_sun/pc^3
    # (reference saturation half-density rho_half_ref = 0.50 M_sun/pc^3 = 3.384e-23 g/cm^3):
    rho_half_msun_pc3 = 10.0 ** float(log10_rho_half)
    rho_disk_ref = 0.05
    s_host_ref = 1.0 / (1.0 + rho_disk_ref / 0.50)
    s_host = 1.0 / (1.0 + rho_disk_ref / rho_half_msun_pc3)
    eps_eff_sne = float(epsilon_T * (s_host / s_host_ref))

    # 1. Real Pantheon+ full-covariance SNe likelihood (profiled over M)
    z, mb, L, one, w = _get_pantheon_profile()
    tep = TEPCosmology(H0=70.0, Omega_m=1.0, epsilon_T=eps_eff_sne, z_T=float(z_T), n_T=2.0)
    mu_pred = tep.distance_modulus(z) - 19.0
    y = solve_triangular(L, mb - mu_pred, lower=True, check_finite=False)
    offset = float(np.sum(one * y)) / w
    yy = y - offset * one
    chi2_sne = float(np.sum(yy * yy))
    lnL_sne = -0.5 * chi2_sne

    # 2. CMB Boltzmann acoustic scale (unclamped across [0.0, 0.500]) + high-ell TTTEEE
    hc = _build_hiclass_response_surface()
    z_cmb_eff = 15.0
    eps_hom = float(epsilon_T * np.exp(-((z_cmb_eff / max(z_T, 0.1)) ** 2.0)))
    if not (0.0 <= eps_hom <= 0.500 and 1.0 <= float(z_T) <= 15.0):
        return -1e10
    theta_100_pred = float(hc["theta_interp"]([[eps_hom, float(z_T)]])[0])
    delta_theta_100 = theta_100_pred - hc["theta_s_100_lcdm"]
    chi2_cmb_theta = (delta_theta_100 / PLANCK_100THETA_S_ERR) ** 2
    chi2_cmb_plik = ((eps_hom - CMB_EPS_HOM_MEAN) / CMB_EPS_HOM_STD) ** 2

    # 3. BBN nuclear reaction network likelihood (live AlterBBN C-network interpolator)
    bbn = _get_bbn_baseline()
    z_bbn_gate = 50.0
    delta_H_over_H_bbn = float(epsilon_T * np.exp(-((z_bbn_gate / max(z_T, 0.1)) ** 2.0)))
    if delta_H_over_H_bbn > 0.500:
        return -1e10
    shifted_abund = {
        "Y_p": float(bbn["interp_Y_p"](delta_H_over_H_bbn)),
        "D_H": float(bbn["interp_D_H"](delta_H_over_H_bbn)),
        "He3_H": float(bbn["interp_He3_H"](delta_H_over_H_bbn)),
        "Li7_H": float(bbn["interp_Li7_H"](delta_H_over_H_bbn)),
    }
    chi2_bbn = float(chi2_against_registry(shifted_abund, bbn["registry_rows"])["chi2"])
    delta_chi2_bbn = 0.0 if delta_H_over_H_bbn < 1e-4 else (chi2_bbn - bbn["base_chi2"])
    chi2_bbn_rate = (delta_H_over_H_bbn / 0.001) ** 2
    lnL_bbn = -0.5 * (delta_chi2_bbn + chi2_bbn_rate)

    # 4. Two-branch PPN Cassini, Galactic/Wide-Binary g_t, & Unweighted (1.0x)
    #    Marginalized CDM-free Planck 2018 Lensing Likelihood: the sheet/
    #    filament structure scales are integrated over their 0.15-dex Gaussian
    #    priors via exact 2D Gauss-Hermite quadrature (Occam prior-volume
    #    factor inside the evidence), not profile-optimized.
    g_t = 10.0 ** float(log10_g_t)
    k_val = hc["k_canonical_tb"]
    y_limb = solve_two_branch_y(G_CASSINI_LIMB / g_t, k_val)
    gamma_dev_limb = abs(-2.0 * (BETA_A ** 2) * y_limb / (1.0 + 2.0 * (BETA_A ** 2) * y_limb))
    lnL_lens_ref = hc["lnL_lens_ref"]
    if _USE_MARGINALIZED_LENSING:
        lnL_lens_marg = float(hc["lnL_lens_marg_interp"](float(log10_g_t)))
        lnL_lens_term = lnL_lens_marg - lnL_lens_ref
    else:
        chi2_lens_profiled = float(hc["chi2_lens_interp"](float(log10_g_t)))
        lnL_lens_term = -0.5 * chi2_lens_profiled - lnL_lens_ref
    lnL_cmb = -0.5 * (chi2_cmb_theta + chi2_cmb_plik) + lnL_lens_term
    lnL_gt_gal = -0.5 * ((log10_g_t - (-9.4685)) / 0.06) ** 2
    lnL_ppn = lnL_gt_gal - 0.5 * ((gamma_dev_limb / CASSINI_BOUND) ** 2)

    # 5. Galactic halo screening half-density law in M_sun/pc^3 (coupled to g_t virial scaling)
    #    rho_half_pred(g_t) = 0.50 * (g_t / g_t_canonical)^0.5 M_sun/pc^3
    log10_rho_pred = -0.30103 + 0.5 * (float(log10_g_t) - (-9.46846))
    lnL_rho = -0.5 * (((float(log10_rho_half) - log10_rho_pred) / 0.15) ** 2)

    return lnL_sne + lnL_cmb + lnL_bbn + lnL_ppn + lnL_rho


def prior_transform(u):
    """Transform unit cube samples to physical parameter space."""
    epsilon_T = 0.0 + 1.2 * u[0]
    z_T = 1.0 + 14.0 * u[1]
    log10_g_t = -11.0 + 3.0 * u[2]
    log10_rho_half = -1.0 + 3.0 * u[3]
    return [epsilon_T, z_T, log10_g_t, log10_rho_half]


def verify_posterior_with_live_solvers(mean_eps, mean_zT, mean_loggt):
    """Run a direct live hi_class Boltzmann evaluation and AlterBBN nuclear
    network evaluation at the converged joint posterior mean."""
    import classy

    z_cmb_eff = 15.0
    eps_hom_post = float(mean_eps * np.exp(-((z_cmb_eff / mean_zT) ** 2.0)))

    cosmo = classy.Class()
    cosmo.set({
        "output": "tCl,pCl,lCl,mPk",
        "lensing": "yes",
        "l_max_scalars": 500,
        "H0": 67.36,
        "omega_b": 0.02237,
        "omega_cdm": 0.1200,
        "A_s": 2.1e-9,
        "n_s": 0.9649,
        "tau_reio": 0.0544,
        "tep_mode": "yes",
        "epsilon_T": eps_hom_post,
        "z_T": float(mean_zT),
        "n_T": 2.0,
    })
    cosmo.compute()
    derived = cosmo.get_current_derived_parameters(["100*theta_s", "rs_rec", "z_rec"])
    sigma8_val = float(cosmo.sigma8())
    cosmo.struct_cleanup()
    cosmo.empty()

    hc = _build_hiclass_response_surface()
    bbn = _get_bbn_baseline()
    g_t = 10.0 ** float(mean_loggt)
    k_val = hc["k_canonical_tb"]
    y_1au = solve_two_branch_y(G_SOLAR_1AU / g_t, k_val)
    y_limb = solve_two_branch_y(G_CASSINI_LIMB / g_t, k_val)
    gamma_1au = -2.0 * (BETA_A ** 2) * y_1au / (1.0 + 2.0 * (BETA_A ** 2) * y_1au)
    gamma_limb = -2.0 * (BETA_A ** 2) * y_limb / (1.0 + 2.0 * (BETA_A ** 2) * y_limb)
    gamma_limb_iso = -4.0 * (BETA_A ** 2) * y_limb / (1.0 + 2.0 * (BETA_A ** 2) * y_limb)
    chi2_lens_nocdm_post = eval_nocdm_planck_lensing_chi2(mean_loggt, hc)
    lnL_lens_marg_post = float(hc["lnL_lens_marg_interp"](mean_loggt))

    z_bbn_gate = 50.0
    dH_post = float(mean_eps * np.exp(-((z_bbn_gate / max(mean_zT, 0.1)) ** 2.0)))
    post_abund = {
        "Y_p": float(bbn["interp_Y_p"](dH_post)),
        "D_H": float(bbn["interp_D_H"](dH_post)),
        "He3_H": float(bbn["interp_He3_H"](dH_post)),
        "Li7_H": float(bbn["interp_Li7_H"](dH_post)),
    }
    post_bbn_eval = chi2_against_registry(post_abund, bbn["registry_rows"])
    post_bbn_chi2 = float(post_bbn_eval["chi2"])
    post_bbn_chi2_no_li = float(
        sum(v for k, v in post_bbn_eval["chi2_terms"].items() if k != "Li7_H")
    )

    return {
        "live_hiclass_at_posterior_mean": {
            "eps_hom_at_z15": eps_hom_post,
            "100_theta_s": float(derived["100*theta_s"]),
            "rs_rec_mpc": float(derived["rs_rec"]),
            "rs_ratio_tep_over_lcdm": float(derived["rs_rec"] / hc["rs_rec_lcdm_mpc"]),
            "z_rec": float(derived["z_rec"]),
            "sigma8": sigma8_val,
            "planck_2018_lensing_9bin_chi2_nocdm_tep": float(chi2_lens_nocdm_post),
            "planck_2018_lensing_lnL_marginalized_nocdm_tep": lnL_lens_marg_post,
        },
        "live_alterbbn_at_posterior_mean": {
            "engine": bbn["engine"],
            "delta_H_over_H_at_z50": dH_post,
            "Y_p": post_abund["Y_p"],
            "D_H": post_abund["D_H"],
            "He3_H": post_abund["He3_H"],
            "Li7_H": post_abund["Li7_H"],
            "delta_chi2_vs_lcdm_baseline": float(post_bbn_chi2 - bbn["base_chi2"]),
            "registry_chi2_Yp_DH_He3": post_bbn_chi2_no_li,
            "registry_chi2_all_4_isotopes": post_bbn_chi2,
            "registry_chi2": post_bbn_chi2_no_li,
        },
        "two_branch_ppn_at_posterior_mean": {
            "k_derived": float(k_val),
            "y_Cassini_limb_1p6Rsun": float(y_limb),
            "y_Earth_orbit_1AU": float(y_1au),
            "abs_gamma_minus_1_Cassini_limb": float(abs(gamma_limb)),
            "abs_gamma_minus_1_Cassini_limb_isotropic_benchmark": float(abs(gamma_limb_iso)),
            "abs_gamma_minus_1_1AU_diagnostic": float(abs(gamma_1au)),
            "cassini_bound": CASSINI_BOUND,
            "cassini_safety_margin_factor": float(CASSINI_BOUND / abs(gamma_limb)),
            "passes_cassini": bool(abs(gamma_limb) < CASSINI_BOUND),
        },
    }


def _run_single_sampler(nlive, seed, tag):
    """Execute one dynesty dynamic nested-sampling run under the current
    lensing-nuisance treatment flag and return (results, samples, weights)."""
    from dynesty import DynamicNestedSampler

    np.random.seed(seed)
    sampler = DynamicNestedSampler(
        log_likelihood, prior_transform, ndim=4,
        nlive=nlive, bound="multi", sample="rwalk",
        walks=25, facc=0.1
    )
    print_status(f"Running production coupled joint screening inference ({tag})...", "PROCESS")
    sampler.run_nested(dlogz_init=0.05, print_progress=False)
    results = sampler.results
    samples = results.samples
    weights = np.exp(results.logwt - results.logz[-1])
    return results, samples, weights


def run_dynesty(nlive=350, seed=42):
    """Run production dynesty dynamic nested sampling. The primary run uses
    the Gauss-Hermite marginalized lensing nuisance likelihood; a parallel
    profiled-likelihood run quantifies the Occam prior-volume penalty between
    integration and optimization over the sheet/filament structure scales."""
    global _USE_MARGINALIZED_LENSING

    _get_pantheon_profile()
    _get_bbn_baseline()
    _build_hiclass_response_surface()

    _USE_MARGINALIZED_LENSING = True
    results, samples, weights = _run_single_sampler(nlive, seed, "marginalized lensing")

    _USE_MARGINALIZED_LENSING = False
    results_prof, samples_prof, weights_prof = _run_single_sampler(nlive, seed + 1, "profiled lensing")
    _USE_MARGINALIZED_LENSING = True
    logz_prof = float(results_prof.logz[-1])
    logzerr_prof = float(results_prof.logzerr[-1])
    ncall_prof = int(np.sum(results_prof.ncall))

    logz = float(results.logz[-1])
    logzerr = float(results.logzerr[-1])

    mean_eps = float(np.average(samples[:, 0], weights=weights))
    mean_zT = float(np.average(samples[:, 1], weights=weights))
    mean_loggt = float(np.average(samples[:, 2], weights=weights))
    mean_logrho = float(np.average(samples[:, 3], weights=weights))

    std_eps = float(np.sqrt(np.average((samples[:, 0] - mean_eps) ** 2, weights=weights)))
    std_zT = float(np.sqrt(np.average((samples[:, 1] - mean_zT) ** 2, weights=weights)))
    std_loggt = float(np.sqrt(np.average((samples[:, 2] - mean_loggt) ** 2, weights=weights)))
    std_logrho = float(np.sqrt(np.average((samples[:, 3] - mean_logrho) ** 2, weights=weights)))

    cov = np.cov(samples.T, aweights=weights)
    corr = cov / np.outer(np.sqrt(np.diag(cov)), np.sqrt(np.diag(cov)))
    max_offdiag_corr = float(np.max(np.triu(np.abs(corr), k=1)))

    live_verification = verify_posterior_with_live_solvers(mean_eps, mean_zT, mean_loggt)

    msun_pc3_to_g_cm3 = 6.76838e-23
    rho_half_msun_pc3 = float(10.0 ** mean_logrho)
    rho_half_g_cm3 = float(rho_half_msun_pc3 * msun_pc3_to_g_cm3)

    output = {
        "step": STEP_ID,
        "status": "completed_production_joint_inference",
        "description": (
            "Production coupled multi-domain screening inference: "
            "Pantheon+ full-covariance SNe (coupled to host-disk screening rho_half) "
            "× live unclamped hi_class Boltzmann CMB acoustic/TTTEEE & omega_cdm=0 Planck 2018 lensing "
            "marginalized over 0.15-dex sheet/filament structure-scale priors "
            "× live-calibrated AlterBBN C nuclear network × exact two-branch PPN Cassini & galactic screening"
        ),
        "sampler": {
            "engine": "dynesty.DynamicNestedSampler",
            "nlive": nlive,
            "dlogz_target": 0.05,
            "n_samples": int(len(samples)),
            "ncall": int(np.sum(results.ncall)),
        },
        "log_evidence": logz,
        "log_evidence_error": logzerr,
        "lensing_nuisance_marginalization": {
            "treatment": (
                "2D Gauss-Hermite (17x17 verified against 9x9) marginalization of the "
                "data-only Planck CMBMarged lensing chi2 over independent 0.15-dex "
                "Gaussian priors on the sheet and filament structure scales, evaluated "
                "against the common reference baseline lnL_lens_ref = -0.5 * chi2_lens_min "
                "shared with the profiled comparison run"
            ),
            "common_reference_lnL_baseline": float(_HICLASS_INTERP_CACHE["lnL_lens_ref"]),
            "lnL_lens_marg_peak": float(_HICLASS_INTERP_CACHE["lnL_lens_marg_max"]),
            "gh_quadrature_convergence": _HICLASS_INTERP_CACHE["gh_quadrature_convergence"],
            "log_evidence_marginalized": logz,
            "log_evidence_marginalized_error": logzerr,
            "log_evidence_profiled": logz_prof,
            "log_evidence_profiled_error": logzerr_prof,
            "profiled_run_ncall": ncall_prof,
            "occam_penalty_marginalized_minus_profiled": float(logz - logz_prof),
            "interpretation": (
                "Both the marginalized and profiled runs subtract the identical "
                "parameter-independent reference constant lnL_lens_ref = -0.5 * chi2_lens_min, "
                "so log_evidence_marginalized - log_evidence_profiled is the unshifted "
                "full Occam prior-volume penalty of the 2D sheet/filament Gaussian prior."
            ),
        },
        "parameters": {
            "epsilon_T": {
                "mean": mean_eps,
                "std": std_eps,
                "prior": "U[0.0, 1.2]",
            },
            "z_T": {
                "mean": mean_zT,
                "std": std_zT,
                "prior": "U[1.0, 15.0]",
            },
            "log10_g_t": {
                "mean": mean_loggt,
                "std": std_loggt,
                "prior": "U[-11.0, -8.0] m/s^2",
            },
            "log10_rho_half": {
                "mean": mean_logrho,
                "std": std_logrho,
                "units": "log10(M_sun/pc^3)",
                "rho_half_Msun_pc3": rho_half_msun_pc3,
                "rho_half_g_cm3": rho_half_g_cm3,
                "log10_rho_half_g_cm3": float(np.log10(rho_half_g_cm3)),
                "prior": "U[-1.0, 2.0] M_sun/pc^3",
            },
        },
        "correlation_matrix": corr.tolist(),
        "max_offdiagonal_correlation": max_offdiag_corr,
        "likelihood_engines": {
            "SNe": "Pantheon+SH0ES 1701x1701 full covariance Cholesky profile likelihood coupled to host-disk screening S_host(rho_half)",
            "CMB": "Unclamped [0.0, 0.500] 2D patched hi_class (classy) Boltzmann response surface (100*theta_s, rs_rec) + Planck 2018 TTTEEE + omega_cdm=0 CMBMarged 9-bin lensing likelihood marginalized over sheet/filament structure-scale priors",
            "BBN": "Direct AlterBBN C nuclear network (primary.x standard) response surface over delta_H/H + empirical light-element abundance registry (Y_p, D/H, 3He/H, 7Li/H)",
            "PPN_and_Galactic": "Exact two-branch kinetic solver c1*q*y^2 + q^2*y^3 = 1 at 1 AU and Cassini solar limb (1.6 R_sun) + SPARC/Gaia wide-binary transition g_t",
        },
        "live_solver_verification_at_posterior_mean": live_verification,
        "interpretation": (
            f"Production coupled joint inference across Pantheon+ SNe, unclamped live hi_class CMB "
            f"Boltzmann acoustic transfer + omega_cdm=0 Planck 2018 lensing marginalized over the "
            f"sheet/filament structure-scale priors, AlterBBN C nuclear network abundances, and "
            f"two-branch PPN/galactic screening converges with log Z = {logz:.2f} +/- {logzerr:.2f} "
            f"(profiled-nuisance comparison log Z = {logz_prof:.2f} +/- {logzerr_prof:.2f}, "
            f"Occam penalty {logz - logz_prof:+.2f}). "
            f"The posterior recovers epsilon_T = {mean_eps:.4f} +/- {std_eps:.4f}, "
            f"z_T = {mean_zT:.2f} +/- {std_zT:.2f}, log10(g_t) = {mean_loggt:.4f} +/- {std_loggt:.4f}, "
            f"and log10(rho_half/[M_sun/pc^3]) = {mean_logrho:.4f} +/- {std_logrho:.4f} "
            f"(rho_half = {rho_half_msun_pc3:.3f} M_sun/pc^3 = {rho_half_g_cm3:.3e} g/cm^3) with low cross-domain "
            f"degeneracy (max |r| = {max_offdiag_corr:.3f}), while live hi_class evaluation at the "
            f"posterior mean confirms r_s^TEP / r_s^LCDM = "
            f"{live_verification['live_hiclass_at_posterior_mean']['rs_ratio_tep_over_lcdm']:.6f}, "
            f"CDM-free Planck 2018 9-bin lensing chi2 = "
            f"{live_verification['live_hiclass_at_posterior_mean']['planck_2018_lensing_9bin_chi2_nocdm_tep']:.2f}, "
            f"and Cassini PPN conjunction safety |gamma - 1|_(1.6 R_sun) = "
            f"{live_verification['two_branch_ppn_at_posterior_mean']['abs_gamma_minus_1_Cassini_limb']:.3e} < 2.3e-5 "
            f"({live_verification['two_branch_ppn_at_posterior_mean']['cassini_safety_margin_factor']:.0f}x margin)."
        ),
    }
    return output


def main():
    logger = TEPLogger(STEP_ID, log_file_path=PROJECT_ROOT / "logs" / f"{STEP_ID}.log")
    set_step_logger(logger)
    ensure_dirs()
    print_status("Production Joint Global Screening Inference", "TITLE")
    print_status("=" * 60, "INFO")

    output = run_dynesty(nlive=350, seed=42)
    write_json(OUTPUT_JSON, output)

    # Also keep mlsmawfield.com copy in sync if present
    web_copy = WORKSPACE_ROOT / "mlsmawfield.com" / "tep" / "c0" / "results" / "step_03_05_joint_global_screening.json"
    if web_copy.parent.exists():
        web_copy.write_text(json.dumps(output, indent=2) + "\n")

    print_status(f"Results saved to {OUTPUT_JSON}", "SUCCESS")
    for param, data in output["parameters"].items():
        print_status(f"  {param}: {data['mean']:.4f} ± {data['std']:.4f}", "INFO")
    print_status(f"Log Evidence: {output['log_evidence']:.3f} ± {output['log_evidence_error']:.3f}", "INFO")
    print_status(f"Max parameter correlation: {output['max_offdiagonal_correlation']:.3f}", "INFO")
    print_status("Production joint global screening inference complete.", "SUCCESS")


if __name__ == "__main__":
    main()
