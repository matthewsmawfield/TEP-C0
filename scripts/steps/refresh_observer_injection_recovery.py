"""Signed noiseless recovery on the current observed-frame likelihood."""
import copy
import json
from pathlib import Path
from refresh_observer_comparison import pipeline

data = pipeline.PantheonData()
assert data.load()
current = json.loads(Path('results/current_model_comparison.json').read_text())
intercept = current['models']['M0a_LCDM']['parameters']['M']
records = {}
for label, amplitude in [('positive', .23), ('negative', -.23)]:
    model = pipeline.ModelTEP(pure_shear=False)
    model.bounds = [tuple(b) for b in model.bounds]
    if amplitude < 0:
        model.bounds[0] = (-1., 1.)
    mock = copy.copy(data)
    mock.mb = model.predict(data.z, {'epsilon_shear_los':amplitude,'M':intercept})+data.observer_correction
    fit, loglike, _ = pipeline.fit_mle(model,mock)
    recovered = float(fit['epsilon_shear_los'])
    error = abs(recovered-amplitude)/abs(amplitude)*100
    assert error < .001, (label, recovered, error)
    records[label] = dict(injected_epsilon_shear_los=amplitude,recovered_epsilon_shear_los=recovered,
                           recovery_error_percent=error,passed=True)
output = dict(scope='Noiseless signed injection recovery of the specified inference model; not external validation of TEP',
              observer_convention='Generated with observer factor and fitted with retained zHEL rows',
              positive_injection_test=records['positive'],negative_injection_test=records['negative'])
Path('results/observer_injection_recovery.json').write_text(json.dumps(output,indent=2)+'\n')
print(json.dumps(output,indent=2))
