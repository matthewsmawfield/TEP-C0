import os
import sys
from pathlib import Path
from cobaya.run import run

# Ensure CLASS is found
project_root = Path("/Users/matthewsmawfield/www/Temporal Equivalence Principle/TEP-C0")
sys.path.insert(0, str(project_root / "external/class/build/lib.macosx-11.1-arm64-cpython-313"))
sys.path.insert(0, str(project_root))
os.environ["TEP_CLASS_PYTHONPATH"] = str(project_root / "external/class/build/lib.macosx-11.1-arm64-cpython-313")

config = {
    "likelihood": {
        "planck_2018_highl_plik.TTTEEE": {"clik_file": str(project_root / "data/external/cobaya_packages/data/planck_2018/baseline/plc_3.0/hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik")},
        "planck_2018_lowl.TT": {"clik_file": str(project_root / "data/external/cobaya_packages/data/planck_2018/baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik")},
        "planck_2018_lowl.EE": {"clik_file": str(project_root / "data/external/cobaya_packages/data/planck_2018/baseline/plc_3.0/low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik")}
    },
    "theory": {"classy": {"extra_args": {"tep_mode": "yes", "lensing": "yes"}}},
    "params": {
        "tep_epsilon_T": {"prior": {"min": 0.0, "max": 2.0}, "ref": 0.89, "proposal": 0.05},
        "tep_z_T": {"prior": {"min": 0.5, "max": 15.0}, "ref": 5.0, "proposal": 0.1},
        "tep_n_T": {"value": 1.0},
        "H0": {"prior": {"min": 30.0, "max": 80.0}, "ref": 50.0, "proposal": 1.0},
        "omega_b": {"prior": {"min": 0.018, "max": 0.025}, "ref": 0.0224, "proposal": 0.0001},
        "omega_cdm": {"value": "lambda H0, omega_b: (H0/100)**2 - omega_b"},
        "tau_reio": {"prior": {"min": 0.01, "max": 0.10}, "ref": 0.054, "proposal": 0.005},
        "logA": {"prior": {"min": 2.5, "max": 3.5}, "ref": 3.044, "proposal": 0.01, "drop": True},
        "A_s": {"value": "lambda logA: 1e-10 * np.exp(logA)"},
        "n_s": {"prior": {"min": 0.90, "max": 1.05}, "ref": 0.966, "proposal": 0.005}
    },
    "sampler": {
        "minimize": {"ignore_prior": False, "max_evals": 1000}
    }
}

try:
    updated_info, sampler = run(config)
    print("MLE found:", sampler.products()["minimum"])
except Exception as e:
    print("Error:", e)
