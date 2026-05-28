import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313')
from classy import Class
c = Class()
c.set({
    'tep_epsilon_T': 0.001, 'tep_z_T': 3.0, 'tep_n_T': 1.0, 
    'H0': 67.5, 'omega_b': 0.0224, 'omega_cdm': 0.12, 'tau_reio': 0.054, 
    'A_s': 2.0989031673191437e-09, 'n_s': 0.966, 'output': 'pCl tCl,pCl,lCl,mPk tCl lCl', 'P_k_max_h/Mpc': 10, 
    'l_max_scalars': 2508, 'lensing': 'yes', 'tep_mode': 'yes', 'non_linear': 'halofit'
})
c.compute()
cls = c.lensed_cl(30)
print('TT Cls up to l=30:', cls['tt'][2])
