import sys
from pathlib import Path
sys.path.insert(0, str(Path("scripts/steps").resolve()))
import numpy as np
from core.tep_cosmology import TEPCosmology
from scripts.steps.step_022_three_model_comparison import PantheonData

data = PantheonData()
data.load()

class ModelTEP3:
    def __init__(self):
        self.param_names = ['H0', 'epsilon_T', 'MB']
        self.n_params = 3
        self.bounds = [(50.0, 100.0), (0.0, 1.0), (-20.5, -18.0)]
        self.z_T = 5.0
    def log_likelihood(self, params, data):
        H0, epsilon_T, MB = params
        Om0 = 1.0  # Fixed to matter-only flat universe
        tep_cosmo = TEPCosmology(H0=H0, Omega_m=Om0, epsilon_T=epsilon_T, z_T=self.z_T)
        mu_tep = tep_cosmo.distance_modulus(data.z)
        residuals = data.mb - (mu_tep + MB)
        try:
            return data.gaussian_loglike(residuals)
        except Exception:
            return -1e10

m1_3 = ModelTEP3()

class neg_logl_class:
    def __init__(self, model, data):
        self.model = model
        self.data = data
    def __call__(self, p):
        return -self.model.log_likelihood(p, self.data)

def custom_fit_mle(model, data):
    from scipy.optimize import minimize
    neg_logl = neg_logl_class(model, data)
    x0 = np.array([67.0, 0.25, -18.1])
    res = minimize(neg_logl, x0, method='L-BFGS-B', bounds=model.bounds, options={'maxiter': 3000})
    return dict(zip(model.param_names, res.x)), -res.fun

m1_3_res, m1_3_nll = custom_fit_mle(m1_3, data)
print("M1 (3 params, Om0=1.0):", m1_3_res, m1_3_nll)

