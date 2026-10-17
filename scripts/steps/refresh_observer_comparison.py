"""Rebuild the matched SN table with zHD distances and zHEL observer factors.

Run from TEP-C0. Keeps the released sample, covariance and existing priors.
Fixed-shape comparisons integrate the bounded linear intercept exactly;
the free-turnover result uses independently refined existing quadrature.
"""
import hashlib
import json
import sys
from pathlib import Path

import dynesty
import numpy as np
from scipy.integrate import quad
from scipy.linalg import solve_triangular, cholesky
from scipy.optimize import minimize, minimize_scalar
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts/steps'))
import step_03_01_three_model_comparison as pipeline


def prepare(data, mask=None):
    if mask is None:
        mask = np.ones(len(data.z), dtype=bool)
    z = data.z[mask]
    mb = (data.mb - data.observer_correction)[mask]
    L = cholesky(data.cov[np.ix_(mask, mask)], lower=True)
    one = solve_triangular(L, np.ones(len(z)), lower=True)
    w = np.sum(one * one)
    norm = -.5 * (2 * np.log(np.diag(L)).sum() + len(z) * np.log(2*np.pi))

    def profile(model, shape):
        params = dict(zip(model.param_names[:-1], np.atleast_1d(shape)))
        params['M'] = -19.
        residual = mb - model.predict(z, params)
        y = solve_triangular(L, residual, lower=True, check_finite=False)
        offset = np.sum(one * y) / w
        y = y - offset * one
        return float(np.sum(y * y)), float(-19. + offset)

    def marginal(model, shape):
        chi, mhat = profile(model, shape)
        prob = ndtr((-16. - mhat)*np.sqrt(w)) - ndtr((-21. - mhat)*np.sqrt(w))
        return float(norm - .5*chi + .5*np.log(2*np.pi/w) - np.log(5.) + np.log(prob))

    return profile, marginal, norm


def scalar_record(model, profile, marginal):
    bounds = model.bounds[0]
    opt = minimize_scalar(lambda x: profile(model, x)[0], bounds=bounds,
                          method='bounded', options={'xatol': 1e-10})
    assert opt.success
    chi, mhat = profile(model, opt.x)
    ref = marginal(model, opt.x)
    value, err = quad(lambda x: np.exp(marginal(model, x)-ref), *bounds,
                      points=[opt.x], epsabs=1e-11, epsrel=1e-7, limit=200)
    return {'chi2': chi, 'logZ': float(ref + np.log(value/(bounds[1]-bounds[0]))),
            'logZ_error': float(err/value), 'parameters': {model.param_names[0]: float(opt.x), 'M': mhat},
            'k': model.n_params, 'method': 'bounded Gaussian intercept integral and adaptive shape quadrature'}


def main():
    data = pipeline.PantheonData()
    assert data.load() and data.redshift_column_used == 'zHD'
    profile, marginal, norm = prepare(data)
    models = {}
    for name, zt in [('M0a_LCDM', None), ('M1_zT5', 5.), ('M1_zT100', 100.)]:
        model = pipeline.ModelLCDM() if zt is None else pipeline.ModelTEP(pure_shear=False)
        if zt is not None:
            model.z_T = zt
        models[name] = scalar_record(model, profile, marginal)
        print(name, models[name], flush=True)
    anchor = json.loads((ROOT/'results/audit_anchored_likelihood_corrected.json').read_text())
    q = json.loads((ROOT/'results/audit_anchored_quadrature.json').read_text())
    assert q['conventions']['observer_redshift_column'] == 'zHEL'
    free = q['runs']['shared']['M1_free_zT']
    models['M1_free_zT'] = {
        'chi2': anchor['M1_free_zT']['shared']['chi2'], 'logZ': free['logZ'],
        'logZ_error': free['last_logZ_change'], 'k': 3,
        'parameters': dict(zip(['epsilon_shear_los', 'z_T', 'M'], anchor['M1_free_zT']['shared']['params'])),
        'method': q['method']}
    models['M2_PureShear'] = dict(models['M0a_LCDM'], method='exact identity with M0a under identical priors')
    for name, model, starts in [
        ('M3_wCDM', pipeline.ModelwCDM(), [[.21, -.70], [.36, -1.]]),
        ('M4_CPL', pipeline.ModelCPL(), [[.45, -.64, -4.99], [.30, -1., 0.]])]:
        bounds = np.array(model.bounds[:-1])
        opts = [minimize(lambda x: profile(model, x)[0], x, bounds=bounds,
                        method='Nelder-Mead', options={'maxiter': 6000, 'xatol': 1e-8, 'fatol': 1e-8}) for x in starts]
        opt = min(opts, key=lambda r: r.fun)
        assert opt.success, opt.message
        chi, mhat = profile(model, opt.x)
        # Same likelihood/prior as sampling M explicitly; integrating M
        # exactly removes an unnecessary numerical sampling dimension.
        sampler = dynesty.NestedSampler(lambda x: marginal(model, x),
                    lambda u: bounds[:, 0] + u*(bounds[:, 1]-bounds[:, 0]),
                    len(bounds), nlive=500, bound='multi', sample='rwalk',
                    rstate=np.random.default_rng(20261011))
        sampler.run_nested(dlogz=.05, print_progress=False)
        models[name] = {'chi2': chi, 'parameters': dict(zip(model.param_names, list(opt.x)+[mhat])),
                        'k': model.n_params, 'logZ': float(sampler.results.logz[-1]),
                        'logZ_error': float(sampler.results.logzerr[-1]),
                        'method': 'exact bounded intercept integral; nested shape sampling, 500 live points, dlogz=0.05, seed=20261011'}
        print(name, models[name], flush=True)
    eds = pipeline.ModelEdS()
    chi, mhat = profile(eds, [])
    models['M0b_EdS'] = {'chi2': chi, 'parameters': {'M': mhat}, 'k': 1,
                          'logZ': marginal(eds, []), 'logZ_error': 0.,
                          'method': 'exact bounded Gaussian intercept integral'}
    for record in models.values():
        record['log_likelihood'] = norm - .5*record['chi2']
        record['BIC'] = record['k']*np.log(len(data.z)) - 2*record['log_likelihood']
        record['BF_vs_LCDM'] = float(np.exp(record['logZ']-models['M0a_LCDM']['logZ']))
        record['delta_chi2_vs_LCDM'] = record['chi2']-models['M0a_LCDM']['chi2']
    source_paths = ['core/cosmology.py', 'scripts/utils/static_metric.py',
                    'scripts/steps/step_03_01_three_model_comparison.py', 'scripts/steps/refresh_observer_comparison.py']
    output = {'sample_size': len(data.z), 'redshift': 'zHD', 'observer_redshift': 'zHEL',
              'luminosity_distance': '(1+zHEL) chi(zHD)', 'exponent': 2, 'H0_reference': 70,
              'M_prior': [-21, -16], 'covariance_sha256': data.provenance['cov_file']['sha256'],
              'models': models, 'source_sha256': {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths}}
    (ROOT/'results/current_model_comparison.json').write_text(json.dumps(output, indent=2)+'\n')
    # Keep low-redshift and calibrator-excision checks on precisely the same convention.
    import pandas as pd
    df = pd.read_csv(data.data_file, sep=r'\s+', comment='#')
    valid = np.isfinite(df.zHD) & np.isfinite(df.m_b_corr) & (df.zHD > .001) & (df.zHD < 3.)
    iscal = df.IS_CALIBRATOR.values[valid].astype(bool)
    cuts = {f'z_above_{cut}': data.z > cut for cut in [.01, .023, .05]}
    cuts['remove_calibrators'] = ~iscal
    checks = {}
    for name, mask in cuts.items():
        p, m, _ = prepare(data, mask)
        lcdm = scalar_record(pipeline.ModelLCDM(), p, m)
        tep = scalar_record(pipeline.ModelTEP(pure_shear=False), p, m)
        checks[name] = {'n': int(mask.sum()), 'LCDM': lcdm, 'TEP': tep,
                        'delta_chi2': tep['chi2']-lcdm['chi2'],
                        'BF_vs_LCDM': float(np.exp(tep['logZ']-lcdm['logZ']))}
        print(name, checks[name]['delta_chi2'], checks[name]['BF_vs_LCDM'], flush=True)
    (ROOT/'results/observer_robustness_quadrature.json').write_text(json.dumps(
        {'conventions': {k:v for k,v in output.items() if k not in ['models', 'source_sha256']}, 'cuts': checks}, indent=2)+'\n')


if __name__ == '__main__':
    main()
