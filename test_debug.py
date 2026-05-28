import sys
sys.path.insert(0, ".")
from cobaya.run import run
from cobaya.theories.classy import classy

original_get_Cl = classy.get_Cl
def patched_get_Cl(self, ell_factor=False, units="FIRASmuK2"):
    cls = original_get_Cl(self, ell_factor=ell_factor, units=units)
    print(f"Cl length: {len(cls['tt'])}")
    return cls
classy.get_Cl = patched_get_Cl

info = "results/outputs/tep_cobaya_sne.updated.yaml"
try:
    run(info, debug=True)
except Exception as e:
    pass
