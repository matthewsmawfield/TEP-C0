from cobaya.model import get_model
import yaml
import numpy as np

with open('results/outputs/tep_cobaya_sne.updated.yaml', 'r') as f:
    info = yaml.safe_load(f)
info["output"] = None
model = get_model(info)

np.random.seed(42)
point = model.prior.sample()
print("Sampled point from priors:")
print({k: float(v) for k, v in point.items()})

try:
    res = model.logposterior(point)
    print("Logpost:", res.logpost)
    for k, v in res.loglikes.items():
        print(f"Likelihood {k}: {v}")
except Exception as e:
    import traceback
    traceback.print_exc()
