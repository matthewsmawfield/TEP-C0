import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313')
from classy import Class
c = Class()
c.set({
    'tep_mode': 'yes',
    'H0': 67.5,
    'omega_b': 0.0224,
    'omega_cdm': 0.120,
    'tau_reio': 0.054,
    'A_s': 2.1e-9,
    'n_s': 0.966,
    'output': 'tCl,pCl,lCl,mPk',
    'P_k_max_h/Mpc': 10,
    'l_max_scalars': 2500,
    'lensing': 'yes',
    'non_linear': 'halofit'
})
c.compute()
print('Halofit Success!')
