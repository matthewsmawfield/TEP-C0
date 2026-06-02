import sys
sys.path.insert(0, "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313")
from classy import Class

def test_no_lambda():
    params = {
        "output": "tCl,pCl,lCl",
        "l_max_scalars": 2500,
        "h": 0.6736,
        "omega_b": 0.02237,
        "omega_cdm": 0.1200,
        "Omega_Lambda": 0.0,  # FORCE Dark Energy to 0
        "A_s": 2.100e-9,
        "n_s": 0.9649,
        "tau_reio": 0.0544,
    }
    cosmo = Class()
    try:
        cosmo.set(params)
        cosmo.compute()
        derived = cosmo.get_current_derived_parameters(["Omega_Lambda"])
        print("SUCCESS!")
        print(derived)
    except Exception as e:
        print("FAILED:", e)

test_no_lambda()
