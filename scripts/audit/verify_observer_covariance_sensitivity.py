"""Matched deterministic checks of the existing covariance controls."""
import copy
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "steps"))
from refresh_observer_comparison import prepare, pipeline, scalar_record

data = pipeline.PantheonData()
assert data.load()
statvar = np.where(data.dmb**2 > 0, data.dmb**2, np.diag(data.cov))
records = {}
for label, covariance in [('full',data.cov), ('diagonal_only',np.diag(np.diag(data.cov))),
                          ('released_error_diagonal',np.diag(statvar)),
                          ('inflated_diagonal_1.5x',np.diag(1.5*np.diag(data.cov)))]:
    variant = copy.copy(data)
    variant.cov = covariance
    profile, marginal, _ = prepare(variant)
    lcdm = scalar_record(pipeline.ModelLCDM(),profile,marginal)
    tep = scalar_record(pipeline.ModelTEP(pure_shear=False),profile,marginal)
    records[label] = dict(LCDM=lcdm,TEP=tep,delta_chi2=tep['chi2']-lcdm['chi2'],
                          BF_vs_LCDM=float(np.exp(tep['logZ']-lcdm['logZ'])))
    print(label,records[label]['delta_chi2'],records[label]['BF_vs_LCDM'],flush=True)
output = dict(convention='zHD distances and zHEL observer factor; unchanged priors',
              released_error_diagonal_scope='m_b_corr_err_DIAG squared with zero-error full-diagonal fallback; not a supplied full statistical covariance',
              variants=records)
Path('results/observer_covariance_sensitivity.json').write_text(json.dumps(output,indent=2)+'\n')
