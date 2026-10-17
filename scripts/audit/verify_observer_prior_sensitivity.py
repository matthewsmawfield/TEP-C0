"""Recompute the manuscript's existing amplitude-prior controls exactly."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "steps"))
from refresh_observer_comparison import prepare, pipeline

data = pipeline.PantheonData()
assert data.load()
profile, marginal, _ = prepare(data)
current = json.loads(Path('results/current_model_comparison.json').read_text())
model = pipeline.ModelTEP(pure_shear=False)
ref = marginal(model, current['models']['M1_zT5']['parameters']['epsilon_shear_los'])
baseline = current['models']['M0a_LCDM']['logZ']
records = []
for low, high in [(0., .5), (0., 1.), (0., 2.), (-1., 1.)]:
    points = [v for v in [.2, .4, .5, .6, .7, .8, .87, .9, 1.] if low < v < high]
    value, error = quad(lambda v: np.exp(marginal(model, v)-ref), low, high,
                        points=points, epsabs=1e-35, epsrel=1e-7, limit=200)
    logz = float(ref + np.log(value/(high-low)))
    records.append(dict(amplitude_prior=[low, high], logZ=logz,
                        BF_vs_LCDM=float(np.exp(logz-baseline)),
                        relative_quadrature_error=float(error/value)))
output = dict(observer_convention='zHD integral with zHEL observer factor',
              fixed_turnover=5., magnitude_prior=[-21., -16.], controls=records)
Path('results/observer_prior_sensitivity.json').write_text(json.dumps(output, indent=2)+'\n')
print(json.dumps(output, indent=2))
