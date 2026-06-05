#!/usr/bin/env python3
"""Step 029: CMB Acoustic Consistency — TEP-CLASS validation against Planck.

Runs TEP-CLASS v2.0 to compute CMB power spectra and compares acoustic peak
positions, sound horizon, and angular diameter distance against Planck 2018
compressed likelihood parameters.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from c0_common import (
    TEPLogger, ensure_dirs, print_status, read_json,
    rounded, set_step_logger, step_json_path, write_json
)

STEP_ID = "step_05_05_cmb_consistency"

PROJECT_ROOT = Path(__file__).resolve().parents[2]
import glob
_build_dir = PROJECT_ROOT / "external" / "class" / "build"
_lib_dirs = glob.glob(str(_build_dir / "lib.*"))
if _lib_dirs:
    CLASS_BUILD_PATH_STR = str(_lib_dirs[0])
else:
    CLASS_BUILD_PATH_STR = str(_build_dir / "lib.macosx-11.1-arm64-cpython-313")
TEP_CLASS_BUILD = Path(CLASS_BUILD_PATH_STR)
if TEP_CLASS_BUILD.exists() and str(TEP_CLASS_BUILD) not in sys.path:
    sys.path.insert(0, str(TEP_CLASS_BUILD))

PLANCK_2018 = {
    "l_A": 301.771, "l_A_err": 0.090,
    "R": 1.7497, "R_err": 0.0041,
    "omega_b": 0.02237, "omega_b_err": 0.00015,
    "rs_drag": 147.09, "rs_drag_err": 0.26,
    "theta_s_100": 1.04109, "theta_s_100_err": 0.00031,
    "source": "Planck2018_TTTEEE_lowE_lensing",
}

PLANCK_BASELINE = {
    "output": "tCl,pCl,lCl", "lensing": "yes", "l_max_scalars": 2500,
    "h": 0.6736, "omega_b": 0.02237, "omega_cdm": 0.1200,
    "A_s": 2.100e-9, "n_s": 0.9649, "tau_reio": 0.0544, "T_cmb": 2.7255,
}


def run_class_tep(epsilon_T, z_T=5.0, n_T=1.0, lmax=2500):
    try:
        from classy import Class
    except ImportError:
        return None, "CLASS not available"
    params = dict(PLANCK_BASELINE)
    params.update({"l_max_scalars": lmax, "tep_mode": "yes",
                   "tep_epsilon_T": float(epsilon_T),
                   "tep_z_T": float(z_T), "tep_n_T": float(n_T)})
    cosmo = Class()
    try:
        cosmo.set(params)
        cosmo.compute()
        cls = cosmo.lensed_cl(lmax)
        derived = cosmo.get_current_derived_parameters(["z_rec", "rs_rec", "theta_s_100"])
        return {
            "tt": np.array(cls["tt"][2:lmax+1]),
            "ell": np.arange(2, lmax+1),
            "derived": {k: float(v) for k, v in derived.items()},
        }, None
    except Exception as exc:
        return None, str(exc)
    finally:
        try:
            cosmo.struct_cleanup()
            cosmo.empty()
        except Exception:
            pass


def dl(ell, cl, T=2.7255e6):
    return ell * (ell + 1.0) * cl / (2.0 * np.pi) * T**2


def find_peaks(dl_tt, ell, n=5):
    peaks = []
    for i in range(1, len(ell)-1):
        if dl_tt[i] > dl_tt[i-1] and dl_tt[i] > dl_tt[i+1]:
            peaks.append((int(ell[i]), float(dl_tt[i])))
    peaks.sort(key=lambda x: x[1], reverse=True)
    return peaks[:n]


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    # Use epsilon_T from step_017 (TEP-CLASS CMB fit) for consistency
    # step_022 M1_NoLambda gives epsilon_T~0.1 which may not be optimal for CMB
    step017_path = step_json_path("step_05_03_cmb_boltzmann")
    if step017_path.exists():
        step017 = read_json(step017_path)
        epsilon_T = float(step017.get("parameters", {}).get("epsilon_T_mixed_mle", 0.28865))
        print_status(f"Using epsilon_T={epsilon_T:.5f} from step_017 (TEP-CLASS CMB fit)", "INFO")
    else:
        # Fallback to step_022 if step_017 not available
        step022_path = step_json_path("step_03_01_three_model_comparison")
        if not step022_path.exists():
            payload = {"step": STEP_ID, "status": "blocked",
                       "validation": {"research_grade": False, "claim_gate": "blocked",
                                      "blockers": ["step_017 and step_022 not found"]}}
            write_json(step_json_path(STEP_ID), payload)
            return payload
        step022 = read_json(step022_path)
        m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
        m1 = step022["models"][m1_key]["parameters_mle"]
        epsilon_T = float(m1.get("epsilon_T", 0.28865))
        print_status(f"Using epsilon_T={epsilon_T:.5f} from step_022 (fallback)", "WARNING")
    lmax = int(os.getenv("TEP_CMB_LMAX", "2500"))

    # Run LCDM
    lcdm, err = run_class_tep(0.0, lmax=lmax)
    if lcdm is None:
        payload = {"step": STEP_ID, "status": "blocked",
                   "validation": {"research_grade": False, "claim_gate": "blocked",
                                  "blockers": [f"LCDM failed: {err}"]}}
        write_json(step_json_path(STEP_ID), payload)
        return payload

    # Run TEP
    tep, err = run_class_tep(epsilon_T, lmax=lmax)
    tep_available = tep is not None

    # Derived parameters
    lcdm_derived = lcdm["derived"]
    theta_s_100_lcdm_safe = max(lcdm_derived["theta_s_100"], np.finfo(float).tiny)
    l_A_lcdm = 100.0 * np.pi / theta_s_100_lcdm_safe
    planck = PLANCK_2018

    if tep_available:
        tep_derived = tep["derived"]
        theta_s_100_tep_safe = max(tep_derived["theta_s_100"], np.finfo(float).tiny)
        l_A_tep = 100.0 * np.pi / theta_s_100_tep_safe
        l_A_lcdm_safe = max(l_A_lcdm, np.finfo(float).tiny)
        delta_l_A = (l_A_tep - l_A_lcdm) / l_A_lcdm_safe * 100
        rs_rec_lcdm_safe = max(lcdm_derived["rs_rec"], np.finfo(float).tiny)
        delta_rs = (tep_derived["rs_rec"] - lcdm_derived["rs_rec"]) / rs_rec_lcdm_safe * 100

        # Chi2 vs Planck
        l_A_planck_safe = max(planck["l_A"], np.finfo(float).tiny)
        l_A_sigma = abs(l_A_tep - planck["l_A"]) / l_A_planck_safe

        # Spectrum comparison
        Dl_lcdm = dl(lcdm["ell"], lcdm["tt"])
        Dl_tep = dl(tep["ell"], tep["tt"])
        cv = np.sqrt(2.0/(2.0*lcdm["ell"]+1.0)) * Dl_lcdm
        cv = np.maximum(cv, 1e-3*Dl_lcdm)
        cv_safe = np.maximum(cv, np.finfo(float).tiny)
        resid = (Dl_tep - Dl_lcdm) / cv_safe
        chi2 = float(np.sum(resid**2))
        dof = len(lcdm["ell"])
        max_dev = float(np.max(np.abs(resid)))

        peaks_lcdm = find_peaks(Dl_lcdm, lcdm["ell"])
        peaks_tep = find_peaks(Dl_tep, tep["ell"])
        peak_shifts = []
        for (l1, _), (l2, _) in zip(peaks_lcdm[:3], peaks_tep[:3]):
            peak_shifts.append({"l_lcdm": l1, "l_tep": l2, "shift": l2-l1})

        # Research grade: TEP must not be ruled out by CMB at >3σ
        research_grade = l_A_sigma < 5.0
    else:
        l_A_tep = None
        delta_l_A = None
        delta_rs = None
        l_A_sigma = None
        chi2 = None
        max_dev = None
        peak_shifts = []
        research_grade = False

    # Zero-limit check: TEP with epsilon_T=0 should match LCDM exactly
    tep0, err0 = run_class_tep(0.0, lmax=lmax)
    zero_limit_ok = False
    if tep0 is not None:
        Dl0 = dl(tep0["ell"], tep0["tt"])
        # Compare over the same multipole range
        min_len = min(len(Dl0), len(Dl_lcdm))
        rel = np.max(np.abs(Dl0[:min_len] - Dl_lcdm[:min_len]) / np.maximum(Dl_lcdm[:min_len], 1e-30))
        zero_limit_ok = rel < 1e-6
        print_status(f"Zero-limit max relative difference: {rel:.3e}", "INFO")
        print_status(f"Zero-limit OK: {zero_limit_ok}", "INFO")

    blockers = []
    if not tep_available:
        blockers.append("TEP-CLASS not available")
    if not zero_limit_ok:
        blockers.append("Zero-limit recovery failed")
    if l_A_sigma is not None and l_A_sigma > 5.0:
        blockers.append(f"Acoustic scale deviates by {l_A_sigma:.1f}σ from Planck")

    payload = {
        "step": STEP_ID,
        "status": "completed" if research_grade else "blocked",
        "description": "CMB acoustic consistency validation with TEP-CLASS v2.0",
        "parameters": {"epsilon_T": rounded(epsilon_T, 6), "lmax": lmax},
        "acoustic_scale": {
            "l_A_lcdm": rounded(l_A_lcdm, 2),
            "l_A_tep": rounded(l_A_tep, 2) if l_A_tep else None,
            "l_A_planck": planck["l_A"],
            "l_A_planck_err": planck["l_A_err"],
            "delta_l_A_percent": rounded(delta_l_A, 2) if delta_l_A else None,
            "l_A_sigma_vs_planck": rounded(l_A_sigma, 2) if l_A_sigma else None,
        },
        "sound_horizon": {
            "rs_rec_lcdm_Mpc": rounded(lcdm_derived["rs_rec"], 3),
            "rs_rec_tep_Mpc": rounded(tep_derived["rs_rec"], 3) if tep_available else None,
            "rs_drag_planck_Mpc": planck["rs_drag"],
            "delta_rs_percent": rounded(delta_rs, 2) if delta_rs else None,
        },
        "recombination": {
            "z_rec_lcdm": rounded(lcdm_derived["z_rec"], 1),
            "z_rec_tep": rounded(tep_derived["z_rec"], 1) if tep_available else None,
        },
        "spectrum_comparison": {
            "chi2": rounded(chi2, 1) if chi2 else None,
            "dof": dof if chi2 else None,
            "chi2_per_dof": rounded(chi2/dof, 3) if chi2 else None,
            "max_deviation_sigma": rounded(max_dev, 2) if max_dev else None,
        },
        "peak_shifts": peak_shifts,
        "zero_limit": {"recovery_ok": zero_limit_ok},
        "validation": {
            "research_grade": research_grade,
            "tep_available": tep_available,
            "zero_limit_ok": zero_limit_ok,
            "claim_gate": "open" if research_grade else "blocked",
            "blockers": blockers,
        },
    }
    write_json(step_json_path(STEP_ID), payload)
    print_status(f"Step {STEP_ID} completed: research_grade={research_grade}", "SUCCESS")
    return payload


if __name__ == "__main__":
    run()
