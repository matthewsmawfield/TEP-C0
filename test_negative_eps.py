import sys
sys.path.insert(0, 'external/class/build/lib.macosx-11.1-arm64-cpython-313')
from classy import Class

cosmo = Class()
cosmo.set({
    'tep_mode': 'yes',
    'tep_epsilon_T': -0.01,
    'tep_z_T': 3.0,
    'output': 'tCl,lCl,mPk'
})
try:
    cosmo.compute()
    print("SUCCESS")
except Exception as e:
    print(f"FAILED: {e}")
