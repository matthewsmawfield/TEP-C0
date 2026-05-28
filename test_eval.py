import sys
sys.path.insert(0, ".")
from cobaya.model import get_model
info = "results/outputs/tep_cobaya_sne.updated.yaml"
model = get_model(info)

# Get the EXACT point that Cobaya tries to evaluate first!
np = __import__('numpy')
np.random.seed(42)

for i in range(5):
    point = {}
    for p, prior in model.prior.priors.items():
        if prior.is_sampled():
            # sample from the prior
            if prior.is_uniform():
                point[p] = np.random.uniform(prior.bounds()[0], prior.bounds()[1])
            else:
                point[p] = prior.reference() # wait, we can just use reference for everything
    
    # Actually, Cobaya gets initial points by drawing from the proposal around the reference!
    point = {}
    for p, prior in model.prior.priors.items():
        if p in model.parameterization.sampled_params():
            point[p] = prior.reference()
            
    print("Evaluating point:")
    print({k: float(v) for k, v in point.items()})
    
    try:
        res = model.logposterior(point)
        print("Logpost:", res.logpost)
        print(res.loglikes)
    except Exception as e:
        print("Error:", e)
    break
