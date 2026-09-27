#!/usr/bin/env python3
"""
TEP-C0 Step 05-11: Acoustic-Gate Branch Consistency
==================================================
The Step 05-10 Jordan-frame diagnostic recovers the CMB acoustic scale
(100*theta_s = 1.0433) on a strictly flat Einstein-de Sitter background
(Omega_m = 1.0, Omega_Lambda = 0.0) at a homogeneous temporal-shear
amplitude epsilon_T^hom = 0.018. The joint Pantheon+ + Planck 2018
TTTEEE MCMC (Step 03-04/03-05) bounds the *same* engine parameter
tep_epsilon_T to epsilon_T^CMB = -0.0015 +/- 0.0037 on the physical,
LambdaCDM-compatible branch.

This step quantifies the branch separation directly from the two recorded
pipeline outputs and classifies the gate honestly: the EdS run is a
counterfactual existence proof that the acoustic scale does not
mathematically require Lambda under the TEP Jordan-frame mapping; it is
not a passed gate on the physical (joint-fit) background, where acoustic
safety is instead supplied by conformal preservation (TEP-HC Boltzmann
closure, r_s^TEP / r_s^LCDM = 0.999994).

Output:
- results/step_05_11_acoustic_gate_consistency.json
"""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

try:
    from scripts.utils.logger import print_status
except Exception:
    def print_status(msg, status="INFO"):
        print(f"[{status}] {msg}")

RESULTS_DIR = PROJECT_ROOT / "results"


def main() -> None:
    jordan_path = RESULTS_DIR / "step_05_10_jordan_frame_proof.json"
    cobaya_path = RESULTS_DIR / "step_03_05_analyze_cobaya.json"
    for p in (jordan_path, cobaya_path):
        if not p.exists():
            raise FileNotFoundError(f"required input missing: {p}")

    jordan = json.loads(jordan_path.read_text())
    cobaya = json.loads(cobaya_path.read_text())

    eps_eds = float(jordan["best_fit"]["epsilon_T"])
    theta_eds = float(jordan["best_fit"]["100_theta_s"])

    eps_cmb = std_cmb = None
    for row in cobaya["constraints"]:
        if row["parameter"] == "tep_epsilon_T":
            eps_cmb = float(row["mean"])
            std_cmb = float(row["std"])
            break
    if eps_cmb is None:
        raise ValueError("tep_epsilon_T constraint not found in Cobaya analysis output")

    tension_sigma = abs(eps_eds - eps_cmb) / std_cmb

    verdict = {
        "step": "step_05_11_acoustic_gate_consistency",
        "description": (
            "Consistency of the EdS acoustic-scale recovery amplitude with the "
            "joint-fit homogeneous shear bound; both are the same tep_epsilon_T "
            "functional of the patched-CLASS engine evaluated on different "
            "background branches"
        ),
        "eds_diagnostic_branch": {
            "source": "step_05_10_jordan_frame_proof.json",
            "background": jordan.get("background"),
            "epsilon_T_hom": eps_eds,
            "100_theta_s": theta_eds,
            "role": "counterfactual existence proof: the acoustic scale does not "
                    "mathematically require Lambda under the TEP Jordan-frame "
                    "mapping",
        },
        "physical_joint_fit_branch": {
            "source": "step_03_05_analyze_cobaya.json",
            "epsilon_T_cmb_mean": eps_cmb,
            "epsilon_T_cmb_std": std_cmb,
            "role": "homogeneous acoustic-sector amplitude measured on the "
                    "LambdaCDM-compatible background selected by "
                    "Pantheon+ + Planck TTTEEE",
        },
        "same_functional": True,
        "functional_note": (
            "Both runs set the identical patched-CLASS parameter tep_epsilon_T; "
            "the amplitudes are the same functional of the same field evaluated "
            "on different background branches and therefore cannot be realized "
            "simultaneously on a single background."
        ),
        "branch_tension_sigma": tension_sigma,
        "gate_classification": (
            "counterfactual_existence_proof_on_eds_branch"
        ),
        "physical_branch_acoustic_safety": (
            "supplied by conformal preservation (TEP-HC hi_class Boltzmann "
            "closure, r_s^TEP / r_s^LCDM = 0.999994), not by the EdS-diagnostic "
            "amplitude"
        ),
        "open_item": (
            "A full Planck TTTEEE confrontation of the no-Lambda homogeneous-"
            "shear branch (not only theta_s) is required to promote the "
            "matter-only acoustic recovery beyond an existence proof; this "
            "depends on the temporal-horizon closure developed in TEP-TH."
        ),
    }

    out_path = RESULTS_DIR / "step_05_11_acoustic_gate_consistency.json"
    out_path.write_text(json.dumps(verdict, indent=2))
    print_status(
        f"epsilon_T^hom (EdS) = {eps_eds}, epsilon_T^CMB = {eps_cmb} +/- {std_cmb}, "
        f"branch separation = {tension_sigma:.2f} sigma"
    )
    print_status(f"wrote {out_path}")


if __name__ == "__main__":
    main()
