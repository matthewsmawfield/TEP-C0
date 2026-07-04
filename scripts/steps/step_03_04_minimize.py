import sys
from pathlib import Path
import os
from typing import Dict, Any

# Ensure correct path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from step_03_04_cobaya_mcmc import create_cobaya_config, check_tep_class_available, check_cobaya_available
from utils.logger import TEPLogger, set_step_logger, print_status

def run_minimize():
    logger = TEPLogger("step_03_04_minimize", log_file_path=Path("logs/step_03_04_minimize.log"))
    set_step_logger(logger)
    print_status("Starting step_03_04_minimize", "TITLE")
    
    output_prefix = str(Path("results/outputs/tep_cobaya_minimize"))
    
    config = create_cobaya_config(
        output_prefix=output_prefix,
        use_planck=True,
        max_samples=10000,
    )
    
    # Replace sampler with minimizer
    config["sampler"] = {
        "minimize": {
            "method": "bobyqa",
            "max_evals": 10000,
            "ignore_prior": False
        }
    }
    
    from cobaya.run import run
    print_status("Running BOBYQA minimizer...", "PROCESS")
    updated_info, sampler = run(config)
    print_status("Minimization complete!", "SUCCESS")
    
if __name__ == "__main__":
    run_minimize()
