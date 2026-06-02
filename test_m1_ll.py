import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import numpy as np
from scripts.steps.step_03_01_three_model_comparison import PantheonData, ModelTEP
data = PantheonData()
data.load()
model = ModelTEP(pure_shear=False, free_z_T=False)
model.z_T = 5.0

print("Evaluating log likelihoods...")
for eps in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
    params = np.array([eps, -19.3])
    ll = model.log_likelihood(params, data)
    print(f"eps={eps:.1f}: logL={ll:.2f}")
