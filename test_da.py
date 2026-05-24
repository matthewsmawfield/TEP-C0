import sys
import numpy as np
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313')
from classy import Class
cosmo = Class()
cosmo.set({'tep_mode': 'yes', 'tep_epsilon_T': 0.1, 'tep_z_T': 3.0, 'output': 'tCl,pCl,lCl', 'lensing': 'yes'})
cosmo.compute()
print("d_a(0.01):", cosmo.angular_distance(0.01))
