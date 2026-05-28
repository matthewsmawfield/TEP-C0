import sys
sys.path.insert(0, ".")
from cobaya.model import get_model
info = {
    "theory": {
        "classy": {
            "path": "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class",
            "extra_args": {"tep_mode": "yes", "non_linear": "halofit"}
        }
    },
    "likelihood": {"planck_2018_lowl.TT": {}},
    "params": {
        "tep_epsilon_T": 0.001, "tep_z_T": 3.0, "tep_n_T": 1.0, 
        "H0": 67.5, "omega_b": 0.0224, "omega_cdm": 0.12, "tau_reio": 0.054, 
        "A_s": 2.0989031673191437e-09, "n_s": 0.966, "A_planck": 1.0
    }
}
try:
    model = get_model(info)
    point = {"tep_epsilon_T": 0.001, "tep_z_T": 3.0, "tep_n_T": 1.0, 
        "H0": 67.5, "omega_b": 0.0224, "omega_cdm": 0.12, "tau_reio": 0.054, 
        "A_s": 2.0989031673191437e-09, "n_s": 0.966, "A_planck": 1.0}
    res = model.logposterior(point)
    print("Commander logpost:", res.logpost)
except Exception as e:
    import traceback
    traceback.print_exc()
