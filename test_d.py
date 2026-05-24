import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313')
import classy
def get_d(eps):
    cosmo = classy.Class()
    cosmo.set({'tep_epsilon_T': eps, 'tep_z_T': 3.0, 'output': 'tCl,lCl,mPk'})
    cosmo.compute()
    return cosmo.luminosity_distance(1.0)
print("eps=0.0:", get_d(0.0))
print("eps=0.1:", get_d(0.1))
