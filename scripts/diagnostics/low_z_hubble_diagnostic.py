#!/usr/bin/env python3
"""
Low-z Hubble-law diagnostic for TEP-C0 Step 03_01.
Prints mu_model - mu_LCDM at key redshifts to audit normalization, sign,
and scale of the temporal shear correction.
"""
import sys, json
from pathlib import Path
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from core.cosmology import CosmologyFLRW
from core.tep_cosmology import TEPCosmology
from core.static_metric import StaticCosmology

RESULTS_JSON = PROJECT_ROOT / "results" / "step_03_01_three_model_comparison.json"

z_test = np.array([0.001, 0.01, 0.05, 0.1, 0.5, 1.0, 2.0])

def main():
    if not RESULTS_JSON.exists():
        print(f"Missing {RESULTS_JSON}")
        return 1

    with open(RESULTS_JSON) as f:
        data = json.load(f)

    models = data.get("models", {})
    m0 = models.get("M0a_LCDM", {}).get("parameters_mle", {})
    m1 = models.get("M1_NoLambda_zT5", {}).get("parameters_mle", {})
    m2 = models.get("M2_PureShear", {}).get("parameters_mle", {})

    H0_ref = 70.0

    # LCDM baseline
    lcdm = CosmologyFLRW(H0=H0_ref, Om0=m0.get("Om0", 0.3))
    mu_lcdm = lcdm.distance_modulus(z_test)

    # TEP mixed model (M1, zT=5)
    tep = TEPCosmology(
        H0=H0_ref, Omega_m=1.0,
        epsilon_T=m1.get("epsilon_T", 0.0),
        z_T=5.0
    )
    mu_tep = tep.distance_modulus(z_test)

    # Pure shear (M2)
    static = StaticCosmology(H0=H0_ref)
    mu_static = static.distance_modulus(z_test)

    print("=" * 80)
    print("Low-z Hubble-law diagnostic")
    print(f"  LCDM:     Om0={m0.get('Om0', 0.3):.4f}")
    print(f"  TEP(M1):  epsilon_T={m1.get('epsilon_T', 0.0):.4f}, z_T=5.0")
    print(f"  Pure(M2): StaticCosmology(H0={H0_ref})")
    print("=" * 80)
    print(f"{'z':>8} {'mu_LCDM':>12} {'mu_TEP':>12} {'mu_Static':>12} {'TEP-LCDM':>12} {'Static-LCDM':>12}")
    print("-" * 80)
    for i, z in enumerate(z_test):
        print(f"{z:8.4f} {mu_lcdm[i]:12.6f} {mu_tep[i]:12.6f} {mu_static[i]:12.6f} {mu_tep[i]-mu_lcdm[i]:12.6f} {mu_static[i]-mu_lcdm[i]:12.6f}")

    # Local Hubble law check: d_L ≈ c*z/H0 at z=0.001
    c = 299792.458  # km/s
    d_l_lcdm = lcdm.luminosity_distance(z_test[0])
    d_l_tep = tep.luminosity_distance(z_test[0])
    d_l_static = static.luminosity_distance(z_test[0])
    d_l_linear = c * z_test[0] / H0_ref
    print()
    print(f"Local Hubble law at z={z_test[0]}")
    print(f"  Linear expectation d_L = c*z/H0 = {d_l_linear:.6f} Mpc")
    print(f"  LCDM    d_L = {d_l_lcdm:.6f} Mpc")
    print(f"  TEP     d_L = {d_l_tep:.6f} Mpc")
    print(f"  Static  d_L = {d_l_static:.6f} Mpc")

    # High-z behaviour: TEP should give larger distances if epsilon_T > 0
    print()
    print(f"High-z behaviour at z={z_test[-1]}")
    print(f"  TEP - LCDM = {mu_tep[-1] - mu_lcdm[-1]:+.6f} mag")
    print(f"  Static - LCDM = {mu_static[-1] - mu_lcdm[-1]:+.6f} mag")

    return 0

if __name__ == "__main__":
    sys.exit(main())
