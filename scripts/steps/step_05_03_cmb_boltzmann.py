#!/usr/bin/env python3
"""Step 017: CMB CLASS resolver.

Runs a vanilla CLASS LCDM reference and then attempts the real TEP CLASS path.
The CMB gate only opens when the active CLASS build accepts TEP background
parameters and validates its Sigma_0=0 LCDM limit.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from c0_common import (
    TEPLogger,
    ensure_dirs,
    print_status,
    read_json,
    rounded,
    set_step_logger,
    step_json_path,
    write_csv,
    write_json,
)
from scripts.utils.cmb_class_resolver import resolve_cmb

STEP_ID = "step_05_03_cmb_boltzmann"


def load_step022() -> dict:
    path = step_json_path("step_03_01_three_model_comparison")
    if not path.exists():
        raise FileNotFoundError("step_022 required")
    return read_json(path)


def load_joint_fit_constraint(parameter: str) -> dict:
    path = step_json_path("step_03_05_analyze_cobaya")
    if not path.exists():
        raise FileNotFoundError("step_03_05_analyze_cobaya required")
    joint = read_json(path)
    for constraint in joint.get("constraints", []):
        if constraint.get("parameter") == parameter:
            return constraint
    raise KeyError(f"joint-fit constraint {parameter} not found")


def run() -> dict:
    logger = TEPLogger(STEP_ID, log_file_path=Path(f"logs/{STEP_ID}.log"))
    set_step_logger(logger)
    print_status(f"Starting {STEP_ID}", "TITLE")
    ensure_dirs()

    step022 = load_step022()
    m1_key = "M1_free_zT" if "M1_free_zT" in step022.get("models", {}) else "M1_NoLambda_zT5"
    m1 = step022["models"][m1_key]["parameters_mle"]
    lmax = int(os.getenv("TEP_CMB_LMAX", "2500"))
    # Physical-branch TEP-CLASS v2.1 parameters from the joint Pantheon+ + Planck
    # fit (step_03_05_analyze_cobaya). The line-of-sight SNe amplitude
    # epsilon_shear_los of the M1 models is a different functional and must not be
    # substituted here (see 9a_parameter_dictionary and step_05_11).
    eps_c = load_joint_fit_constraint("tep_epsilon_T")
    zt_c = load_joint_fit_constraint("tep_z_T")
    epsilon_t = float(eps_c["mean"])
    z_t = float(zt_c["mean"])
    cmb_run = resolve_cmb(
        epsilon_t=epsilon_t,
        z_t=z_t,
        n_t=1.0,  # Power-law index, fixed in the joint chain
        lmax=lmax,
    )

    if cmb_run.lcdm_reference:
        write_csv(Path(f"results/outputs/{STEP_ID}_class_lcdm_samples.csv"), cmb_run.lcdm_reference["samples"])
        write_csv(Path(f"results/outputs/{STEP_ID}_class_lcdm_tt.csv"), cmb_run.lcdm_reference["tt_power"])
        print_status("CLASS LCDM reference spectrum completed", "SUCCESS")
    if cmb_run.tep_spectra:
        write_csv(Path(f"results/outputs/{STEP_ID}_tep_tt.csv"), cmb_run.tep_spectra["tt_power"])
        print_status("TEP CLASS spectrum completed", "SUCCESS")
    else:
        print_status("TEP CLASS spectrum unavailable in active classy build", "WARNING")

    payload = {
        "step": STEP_ID,
        "status": "completed" if cmb_run.validation.get("research_grade_cmb") else "blocked",
        "description": "CLASS-backed CMB resolver with strict TEP parameter detection",
        "tier": "Resolver",
        "parameters": {
            "H0_mixed_mle": rounded(70.0, 4),  # dimensionless model fixes H0_ref
            "Om0_mixed_mle": rounded(m1.get("Om0", 1.0), 6),  # M1_NoLambda is matter-only
            "Sigma_0_mixed_mle": rounded(m1.get("Sigma_0", 0.0), 8),  # M1 has no Sigma_0
            "tep_epsilon_T": rounded(epsilon_t, 8),
            "tep_epsilon_T_source": "step_03_05_analyze_cobaya joint Pantheon+ + Planck TTTEEE posterior mean",
            "tep_z_T": rounded(z_t, 4),
            "tep_z_T_source": "step_03_05_analyze_cobaya joint posterior mean (unconstrained; spectrum insensitive at epsilon_T ~ 0)",
            "tep_n_T": 1.0,
            "lmax": lmax,
        },
        "cmb_results": {
            "class_lcdm_reference": cmb_run.lcdm_reference,
            "tep_zero_limit": cmb_run.tep_zero_limit,
            "tep_spectra": cmb_run.tep_spectra,
        },
        "validation": cmb_run.validation,
    }
    write_json(step_json_path(STEP_ID), payload)
    return payload


if __name__ == "__main__":
    result = run()
    print({
        "step": result["step"],
        "tier": result["tier"],
        "validation": result["validation"],
    })
