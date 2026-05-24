import sys
import os
from cobaya.model import get_model
sys.path.insert(0, '/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0')
import scripts.steps.step_033_cobaya_tep_inference as s33
config = s33.create_cobaya_config("results/outputs/test_eval", use_planck=False, max_samples=1)
config["resume"] = False
config["force"] = True
model = get_model(config)
point = {
    "tep_epsilon_T": 0.05,
    "tep_z_T": 3.0,
    "H0": 67.4,
    "omega_b": 0.0224,
    "omega_cdm": 0.12,
    "tau_reio": 0.054,
    "A_s": 2.1e-9,
    "n_s": 0.966
}
logpost, loglikes, _ = model.logposterior(point)
print("Log Posterior:", logpost)
print("Log Likelihoods:", loglikes)
