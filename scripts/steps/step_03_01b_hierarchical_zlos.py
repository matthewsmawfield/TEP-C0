#!/usr/bin/env python3
"""Turnover-prior sensitivity from the current full-covariance likelihood.

Integrates the existing LogUniform(1,200) alternative directly, rather than
interpolating three historical Bayes factors. The magnitude intercept has
its exact bounded Gaussian integral; quadrature refinement measures the
turnover integration accuracy. This alternative does not replace the
uniform-turnover prior used in the primary comparison.
"""
import sys
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'scripts'), str(ROOT/'scripts/steps')]
from refresh_observer_comparison import prepare, scalar_record
from step_03_01_three_model_comparison import PantheonData, ModelLCDM, ModelTEP


def main():
    data = PantheonData()
    assert data.load() and data.redshift_column_used == 'zHD'
    profile, marginal, _ = prepare(data)
    lcdm = scalar_record(ModelLCDM(), profile, marginal)
    model = ModelTEP(pure_shear=False)
    # Center on the supported fixed-5 evidence to keep integrals well scaled.
    ref = scalar_record(model, profile, marginal)['logZ']
    lo, hi = np.log(1.), np.log(200.)
    levels = []
    for order in (32, 64, 128):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        terms = []
        for t, weight in zip(lo + (nodes+1)*(hi-lo)/2, weights):
            model.z_T = float(np.exp(t))
            val, err = quad(lambda eps: np.exp(marginal(model, eps)-ref), 0., 2.,
                            points=[.2,.4,.6,.8,1.,1.2], epsabs=1e-10, epsrel=2e-7)
            terms.append(weight*val/2.)
        logz = float(ref+np.log(sum(terms)/2.))
        levels.append({'order': order, 'logZ': logz,
                       'BF_vs_LCDM': float(np.exp(logz-lcdm['logZ']))})
        print(levels[-1], flush=True)
    change = abs(levels[-1]['logZ']-levels[-2]['logZ'])
    assert change < 1e-5, 'Turnover quadrature has not converged'
    output = {'conventions': {'distance_redshift': 'zHD', 'observer_redshift': 'zHEL',
                             'n_sne': len(data.z), 'M_prior': [-21.,-16.], 'epsilon_prior': [0.,2.]},
              'hyperprior': 'LogUniform(1, 200)', 'method': 'exact bounded intercept integral; adaptive amplitude and refined log-turnover quadrature',
              'refinements': levels, 'last_logZ_change': change,
              'BF_marginalized': levels[-1]['BF_vs_LCDM'],
              'lnBF_marginalized': levels[-1]['logZ']-lcdm['logZ'],
              'interpretation': 'Prior sensitivity alternative; the primary free-turnover comparison retains Uniform(0.1,150).'}
    (ROOT/'results/step_03_01b_hierarchical_zlos.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__ == '__main__':
    main()
