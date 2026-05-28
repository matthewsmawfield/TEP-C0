import sys
import yaml
import numpy as np
from cobaya.model import get_model

info = {
    "params": {
        "tep_epsilon_T": 0.001,
        "tep_z_T": 3.0,
        "tep_n_T": 1.0,
        "H0": 67.5,
        "omega_b": 0.0224,
        "omega_cdm": 0.120,
        "tau_reio": 0.054,
        "logA": 3.044,
        "A_s": {"value": "lambda logA: 1e-10 * __import__('numpy').exp(logA)"},
        "n_s": 0.966,
    },
    "likelihood": {
        "core.pantheon_cobaya_likelihood.PantheonCobaya": {}
    },
    "theory": {
        "classy": {
            "path": "/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0/external/class",
            "extra_args": {
                "N_ur": 2.0328, "N_ncdm": 1, "m_ncdm": 0.06,
                "output": "tCl,pCl,lCl,mPk", "P_k_max_h/Mpc": 10,
                "l_max_scalars": 2500, "lensing": "yes", "tep_mode": "yes"
            }
        }
    }
}
info["output"] = None
model = get_model(info)
print("Testing Pantheon only...")
point = {
    "tep_epsilon_T": 0.001, "tep_z_T": 3.0, "tep_n_T": 1.0,
    "H0": 67.5, "omega_b": 0.0224, "omega_cdm": 0.120,
    "tau_reio": 0.054, "logA": 3.044, "n_s": 0.966
}
res = model.logposterior(point)
print("Logpost:", res.logpost)
