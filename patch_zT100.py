import sys
import json
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "scripts" / "steps"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from step_03_01_three_model_comparison import *

def run_patch():
    data = PantheonData()
    data.load()

    model = ModelTEP(pure_shear=False)
    model.z_T = 100.0
    model_name = "M1_Unscreened_zT100"

    print("Fitting MLE...")
    mle, logl, _ = fit_mle(model, data)
    residuals = data.mb - model.predict(data.z, mle, data)
    chi2 = data.chi2(residuals)
    n = len(data.z)
    k = model.n_params
    deviance = chi2 + data.cov_logdet + n * np.log(2.0 * np.pi)
    aic = deviance + 2*k
    bic = deviance + k * np.log(n)

    payload = {
        'parameters_mle': mle,
        'log_likelihood_mle': float(logl),
        'chi2_mle': float(chi2),
        'chi2_red_mle': float(chi2 / (n - k)),
        'aic': float(aic),
        'bic': float(bic),
        'n_params': int(k),
    }

    print(f"MLE: {mle}, logL: {logl}")
    print("Running nested sampling...")
    nested = run_nested_evidence(model, data, nlive=500, dlogz=0.1)
    payload.update(nested)
    
    print("Updating JSON...")
    json_path = Path('results/step_03_01_three_model_comparison.json')
    with open(json_path, 'r') as f:
        res = json.load(f)

    res['models'][model_name] = payload

    m0_ev = res['models']['M0a_LCDM']['log_evidence']
    ln_bf = payload['log_evidence'] - m0_ev
    res['bayes_factors'][f'ln_BF_{model_name}_vs_M0a_LCDM'] = float(ln_bf)
    res['bayes_factors'][f'BF_{model_name}_vs_M0a_LCDM'] = float(np.exp(ln_bf))
    
    if ln_bf > 5:
        interp = "Decisive"
    elif ln_bf > 2.5:
        interp = "Strong"
    elif ln_bf > 1:
        interp = "Substantial"
    else:
        interp = "Weak"
    res['bayes_factors'][f'interpretation_{model_name}_vs_M0a_LCDM'] = interp

    with open(json_path, 'w') as f:
        json.dump(res, f, indent=2)

    print(f"Done. ln BF = {ln_bf:.2f} ({interp})")

if __name__ == "__main__":
    run_patch()
