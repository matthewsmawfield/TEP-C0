import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313")
from classy import Class

def run_test(epsilon_T, z_T=5.0):
    params = {
        "output": "tCl,pCl,lCl",
        "l_max_scalars": 2500,
        "lensing": "yes",
        "h": 0.6736,
        "omega_b": 0.02237,
        # Force Omega_m = 1 by setting omega_cdm = h^2 - omega_b
        "omega_cdm": 0.6736**2 - 0.02237, 
        "A_s": 2.100e-9,
        "n_s": 0.9649,
        "tau_reio": 0.0544,
        "tep_mode": "yes",
        "tep_epsilon_T": epsilon_T,
        "tep_z_T": z_T,
        "tep_n_T": 1.0,
    }
    
    cosmo = Class()
    try:
        cosmo.set(params)
        cosmo.compute()
        cls = cosmo.lensed_cl(2500)
        derived = cosmo.get_current_derived_parameters(["z_rec", "rs_rec", "theta_s_100", "Omega_Lambda"])
        return derived, cls
    except Exception as e:
        return str(e), None

print("Running EdS (Omega_m=1, eps=0)...")
res0 = run_test(0.0)
print("Omega_Lambda:", res0[0].get("Omega_Lambda"))
print("theta_s_100:", res0[0].get("theta_s_100"))

print("\nRunning TEP (Omega_m=1, eps=0.89)...")
resTEP = run_test(0.89)
if isinstance(resTEP[0], dict):
    print("Omega_Lambda:", resTEP[0].get("Omega_Lambda"))
    print("theta_s_100:", resTEP[0].get("theta_s_100"))
else:
    print("Error:", resTEP[0])
