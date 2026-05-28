import sys
sys.path.insert(0, ".")
from cobaya.model import get_model
from cobaya.theories.classy import classy

# Monkey patch!
_old_get_can_support_params = classy.get_can_support_params
def _new_get_can_support_params(self):
    return _old_get_can_support_params(self) + ["tep_epsilon_T", "tep_z_T", "tep_n_T"]
classy.get_can_support_params = _new_get_can_support_params

import yaml
with open('results/outputs/tep_cobaya_sne.updated.yaml', 'r') as f:
    info = yaml.safe_load(f)

# Force a test point
info["output"] = None
model = get_model(info)
print("Model built!")

# Does it pass tep_epsilon_T to classy?
# We can check by looking at the parameters that the classy provider receives!
params = model.parameterization.prior.reference()
print("Evaluating point...")
try:
    res = model.logposterior(params)
    print("Logpost:", res.logpost)
except Exception as e:
    import traceback
    traceback.print_exc()
