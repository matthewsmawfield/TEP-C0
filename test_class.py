import sys
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class/build/lib.macosx-11.1-arm64-cpython-313')
import classy
print("Classy loaded")
cosmo = classy.Class()
cosmo.set({'tep_epsilon_T': 0.1, 'tep_z_T': 2.0, 'output': 'tCl,lCl,mPk'})
cosmo.compute()
z = 1.0
dl = cosmo.luminosity_distance(z)
da = cosmo.angular_distance(z)
print(f"z={z}, dL={dl}, dA={da}, ratio={dl/(da*(1+z)**2)}")
