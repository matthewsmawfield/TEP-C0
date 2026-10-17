"""Deterministic evidence check, with the linear intercept integrated exactly.

Uses the corrected anchored audit's controlled sample and uniform priors.
Gauss-Legendre refinement checks the remaining turnover integral. No Monte
Carlo noise or extra Cepheid error diagonal is introduced.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import log_ndtr
from scipy.linalg import solve_triangular

# Load definitions without executing the expensive nested-sampling runner.
source = Path('scripts/audit/audit_anchored_likelihood.py').read_text()
exec(compile(source.split('# --- MLE comparison')[0], 'audit_definitions', 'exec'))
one = solve_triangular(cholesky(cov, lower=True), np.ones(len(mb)), lower=True)
L = cholesky(cov, lower=True)
w = float(np.sum(one * one))
norm = -0.5*(2*np.log(np.diag(L)).sum()+len(mb)*np.log(2*np.pi))
saved = json.loads(Path('results/audit_anchored_likelihood_corrected.json').read_text())

def marginalized(model, names, shape, residual):
    # Center around M=-19 to avoid subtracting large nearly equal chi squares.
    r = residual(model, names, list(shape)+[-19.])
    y = solve_triangular(L, r, lower=True, check_finite=False)
    offset = float(np.sum(one * y))/w
    yy = y-offset*one
    mhat = -19.+offset
    hi, lo = (-16.-mhat)*np.sqrt(w), (-21.-mhat)*np.sqrt(w)
    logprob = log_ndtr(hi) + np.log1p(-np.exp(log_ndtr(lo)-log_ndtr(hi)))
    return norm - .5*float(np.sum(yy * yy))+.5*np.log(2*np.pi/w)-np.log(5.)+logprob

output = {'conventions': saved['conventions'], 'method': 'exact bounded Gaussian intercept integral; adaptive shape quadrature; refined turnover quadrature', 'runs': {}}
for variant, residual in [('shared', resid_shared), ('anchored', resid_anchored)]:
    logs = {}
    for name, model, names, bounds, x0 in MODELS:
        ref = norm-.5*saved[name][variant]['chi2']
        def shape_integral(zturn=None):
            lo, hi = bounds[0]
            def f(v):
                shape = [v] if zturn is None else [v,zturn]
                return np.exp(marginalized(model,names,shape,residual)-ref)
            points = [x for x in (.2,.3,.4,.6,.7,.8,.9,1.,1.1,1.3) if lo<x<hi]
            value, error = quad(f,lo,hi,points=points,epsabs=1e-10,epsrel=2e-7,limit=200)
            return value/(hi-lo), error/(hi-lo)
        if len(names)==1:
            val, err = shape_integral()
            record = {'logZ': float(ref+np.log(val)), 'relative_quad_error':float(err/val)}
        else:
            levels = []
            # Uniform z_T prior, integrated in log(z_T) with its Jacobian.
            zl, zh = bounds[1]; tl, th = np.log(zl), np.log(zh)
            for order in (32,64,128):
                nodes, weights = np.polynomial.legendre.leggauss(order)
                terms = []
                for t, weight in zip(tl+(nodes+1)*(th-tl)/2, weights):
                    zt = np.exp(t)
                    val, _ = shape_integral(zt)
                    terms.append(weight*zt*val)
                total = sum(terms)*(th-tl)/2/(zh-zl)
                levels.append({'order':order,'logZ':float(ref+np.log(total))})
            record = {'logZ':levels[-1]['logZ'],'turnover_refinements':levels,
                      'last_logZ_change':abs(levels[-1]['logZ']-levels[-2]['logZ'])}
        logs[name]=record
        print(variant,name,json.dumps(record),flush=True)
    baseline=logs['M0a_LCDM']['logZ']
    for name, record in logs.items():
        record['dlnZ_vs_LCDM']=record['logZ']-baseline
        record['BF_vs_LCDM']=float(np.exp(record['dlnZ_vs_LCDM']))
    output['runs'][variant]=logs
    Path('results/audit_anchored_quadrature.json').write_text(json.dumps(output,indent=2)+'\n')
