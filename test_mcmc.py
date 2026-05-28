from cobaya.run import run
import sys
sys.path.insert(0, "core")
from pantheon_cobaya_likelihood import PantheonCobaya

info = {
    "params": {
        "tep_epsilon_T": {"prior": {"min": 0.0, "max": 0.4}, "ref": 0.02, "proposal": 0.005},
        "omega_b": 0.022,
        "omega_cdm": 0.12,
        "H0": 67.0,
    },
    "likelihood": {
        "sn": {"external": PantheonCobaya}
    },
    "theory": {
        "classy": {
            "extra_args": {"tep_mode": "yes"}
        }
    },
    "sampler": {"evaluate": {}}
}

updated_info, sampler = run(info)
print("Evaluate successful!")
