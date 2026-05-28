from cobaya.model import get_model
info = {
    "params": {
        "tep_epsilon_T": 0.2,
        "omega_b": 0.022,
        "omega_cdm": 0.12,
        "H0": 67.0,
    },
    "likelihood": {
        "dummy": {
            "external": lambda _self, tep_epsilon_T, Hubble: -0.5 * (tep_epsilon_T - 0.2)**2,
            "requires": {"tep_epsilon_T": None, "Hubble": {"z": [0.0]}}
        }
    },
    "theory": {
        "classy": {
            "ignore_obsolete": True,
            "extra_args": {
                "tep_mode": "yes"
            }
        }
    }
}
model = get_model(info)
print("Model built successfully!")
model.logposterior({})
print("Logpost computed successfully!")
