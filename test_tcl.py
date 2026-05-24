import sys
import numpy as np
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313')
from classy import Class
cosmo = Class()
cosmo.set({'tep_mode': 'yes', 'tep_epsilon_T': 0.1, 'tep_z_T': 3.0, 'H0': 70.0, 'omega_b': 0.0224, 'omega_cdm': 0.12, 'A_s': 2.1e-9, 'n_s': 0.96, 'output': 'tCl,pCl,lCl', 'lensing': 'yes'})
cosmo.compute()
cl = cosmo.lensed_cl(2500)
for k in cl:
    if k != 'ell':
        has_nan = np.any(np.isnan(cl[k]))
        print(f"{k} has NaN:", has_nan)
