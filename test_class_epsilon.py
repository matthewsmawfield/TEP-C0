import sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path("external/class/build/lib.macosx-11.1-arm64-cpython-313")))
from classy import Class

for eps in [0.0, 0.01, 0.05, -0.05]:
    try:
        cosmo = Class()
        cosmo.set({
            "output": "tCl,pCl,lCl", "l_max_scalars": 2500, "lensing": "yes",
            "h": 0.6736, "omega_b": 0.02237, "omega_cdm": 0.1200,
            "A_s": 2.100e-9, "n_s": 0.9649, "tau_reio": 0.0544,
            "tep_mode": "yes", "tep_epsilon_T": eps, "tep_z_T": 3.0, "tep_n_T": 1.0
        })
        cosmo.compute()
        cls = cosmo.lensed_cl(2500)
        tt = cls["tt"]
        print(f"eps={eps}: TT has NaN: {np.isnan(tt).any()}, TT has inf: {np.isinf(tt).any()}, max_TT: {np.max(tt)}")
    except Exception as e:
        print(f"eps={eps}: Exception {type(e).__name__}: {e}")
