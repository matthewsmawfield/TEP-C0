"""Compare shared-intercept and Cepheid-anchored likelihoods.

Corrected 10 October 2026: use independent CEPH_DIST rather than the
SN-derived MU_SH0ES. Supplied covariance already includes Cepheid-host
uncertainty; do not add the plotting-only MU_SH0ES_ERR_DIAG again.
Retains corpus z>0.001 sample for a controlled comparison, unlike the
official zHD>0.01-or-calibrator sample. Historical output is preserved.
"""
import sys, json, os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.linalg import cholesky, solve_triangular

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts' / 'steps'))
import importlib.util
spec = importlib.util.spec_from_file_location('s301', 'scripts/steps/step_03_01_three_model_comparison.py')
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

NLIVE = int(os.getenv('TEP_ANCHOR_NLIVE', '200'))
DLOGZ = float(os.getenv('TEP_ANCHOR_DLOGZ', '0.3'))

data = m.PantheonData(); assert data.load()
z, mb, cov = data.z, data.mb, data.cov
df = pd.read_csv(data.data_file, sep=r'\s+', comment='#')
valid = np.isfinite(df[data.redshift_column_used].values) & np.isfinite(df['m_b_corr'].values) \
        & (df[data.redshift_column_used].values > 0.001) & (df[data.redshift_column_used].values < 3.0)
isc = df['IS_CALIBRATOR'].values[valid].astype(bool)
cepheid_distance = df['CEPH_DIST'].values[valid]
assert np.isfinite(cepheid_distance[isc]).all()
assert len(isc) == len(mb)

# Release covariance already includes Cepheid-host uncertainty.
cov_a = cov.copy()

def make_gauss(c):
    L = cholesky(c, lower=True, check_finite=False)
    logdet = float(2.0 * np.sum(np.log(np.diag(L))))
    n = len(mb)
    def chi2(r):
        y = solve_triangular(L, r, lower=True, check_finite=False)
        return float(np.sum(y * y))
    def loglike(r):
        return -0.5 * (chi2(r) + logdet + n * np.log(2.0 * np.pi))
    return chi2, loglike

chi2_s, ll_s = make_gauss(cov)      # shared-intercept construction
chi2_a, ll_a = make_gauss(cov_a)    # anchored construction

mlcdm = m.ModelLCDM(); mtep = m.ModelTEP(pure_shear=False)
mtep_free = m.ModelTEP(pure_shear=False, free_z_T=True)

def pred_field(model, params_dict):
    return model.predict(z, {**params_dict, 'M': 0.0}, data) + data.observer_correction

def resid_shared(model, shape_names, x):
    """Residuals under the shared-intercept construction (all rows use model mu + M)."""
    ns = len(shape_names)
    mu = pred_field(model, dict(zip(shape_names, x[:ns]))) + x[ns]
    return mb - mu

def resid_anchored(model, shape_names, x):
    """Residuals under the anchored construction.

    Calibrator rows: m_b - CEPH_DIST - M (independent of the cosmological model).
    Field rows:      m_b - mu_model(z) - M.
    """
    ns = len(shape_names)
    M = x[ns]
    r = np.empty_like(mb)
    r[isc] = mb[isc] - cepheid_distance[isc] - M
    mu_f = pred_field(model, dict(zip(shape_names, x[:ns])))
    r[~isc] = mb[~isc] - mu_f[~isc] - M
    return r

def fit_mle(model, shape_names, bounds, resid, chi2fun, x0):
    def nll(x):
        try:
            return chi2fun(resid(model, shape_names, x))
        except Exception:
            return 1e10
    starts = [x0]
    if 'z_T' in shape_names:
        starts += [[x0[0], zz, x0[-1]] for zz in (0.5, 2., 5., 20., 100.)]
    fits = [minimize(nll, start, method='L-BFGS-B', bounds=bounds, options={'maxiter': 3000}) for start in starts]
    res = min(fits, key=lambda fit: fit.fun)
    assert res.success, res.message
    return float(res.fun), [float(v) for v in res.x]

MODELS = [
    ('M0a_LCDM',        mlcdm,     ['Om0'],                        [(0.05, 0.9), (-21, -16)], [0.35, -19.34]),
    ('M1_zT5',          mtep,      ['epsilon_shear_los'],          [(0.0, 2.0), (-21, -16)],  [0.87, -19.33]),
    ('M1_free_zT',      mtep_free, ['epsilon_shear_los', 'z_T'],   [(0.0, 2.0), (0.1, 150.0), (-21, -16)], [0.87, 76.0, -19.33]),
]

out = {'conventions': {
    'redshift_column': data.redshift_column_used,
    'observer_redshift_column': 'zHEL',
    'luminosity_distance': '(1+zHEL) chi(z_distance)',
    'n_calibrators': int(isc.sum()),
    'anchored_residual_calibrator': 'm_b - CEPH_DIST - M',
    'anchored_residual_field': 'm_b - mu_model(z) - M',
    'anchored_covariance': 'supplied full Pantheon+SH0ES covariance, without additional diagonal inflation',
    'nested_nlive': NLIVE, 'nested_dlogz': DLOGZ,
}}

# --- MLE comparison under both constructions ---
for name, model, snames, bounds, x0 in MODELS:
    c_s, p_s = fit_mle(model, snames, bounds, resid_shared, chi2_s, x0)
    c_a, p_a = fit_mle(model, snames, bounds, resid_anchored, chi2_a, x0)
    out[name] = {
        'shared':  {'chi2': c_s, 'params': p_s},
        'anchored': {'chi2': c_a, 'params': p_a},
    }

out['mle_comparison'] = {
    'dchi2_TEP_zT5_vs_LCDM_shared':   out['M1_zT5']['shared']['chi2']   - out['M0a_LCDM']['shared']['chi2'],
    'dchi2_TEP_zT5_vs_LCDM_anchored': out['M1_zT5']['anchored']['chi2'] - out['M0a_LCDM']['anchored']['chi2'],
    'dchi2_TEP_free_vs_LCDM_shared':  out['M1_free_zT']['shared']['chi2']  - out['M0a_LCDM']['shared']['chi2'],
    'dchi2_TEP_free_vs_LCDM_anchored': out['M1_free_zT']['anchored']['chi2'] - out['M0a_LCDM']['anchored']['chi2'],
}
print(json.dumps(out['mle_comparison'], indent=2), flush=True)

# --- Nested-sampling evidence under both constructions ---
try:
    if os.getenv('TEP_ANCHOR_MLE_ONLY') == '1':
        raise RuntimeError('Nested sampling explicitly skipped; no new BF claimed')
    import dynesty
    def run_nested(model, snames, bounds, resid, llfun):
        b = np.array(bounds, dtype=float)
        lo, w = b[:, 0], b[:, 1] - b[:, 0]
        def pt(u): return lo + u * w
        def ll(u):
            return llfun(resid(model, snames, u))
        s = dynesty.NestedSampler(ll, pt, len(bounds), nlive=NLIVE)
        s.run_nested(dlogz=DLOGZ, print_progress=False)
        return float(s.results.logz[-1]), float(s.results.logzerr[-1])

    for name, model, snames, bounds, x0 in MODELS:
        lz_s, e_s = run_nested(model, snames, bounds, resid_shared, ll_s)
        lz_a, e_a = run_nested(model, snames, bounds, resid_anchored, ll_a)
        out[name]['shared'].update({'logZ': lz_s, 'logZ_err': e_s})
        out[name]['anchored'].update({'logZ': lz_a, 'logZ_err': e_a})
        print(f"{name}: shared lnZ={lz_s:.2f}+-{e_s:.2f}  anchored lnZ={lz_a:.2f}+-{e_a:.2f}", flush=True)

    out['evidence_comparison'] = {
        'dlnZ_zT5_shared':   out['M1_zT5']['shared']['logZ']   - out['M0a_LCDM']['shared']['logZ'],
        'dlnZ_zT5_anchored': out['M1_zT5']['anchored']['logZ'] - out['M0a_LCDM']['anchored']['logZ'],
        'dlnZ_free_shared':  out['M1_free_zT']['shared']['logZ']  - out['M0a_LCDM']['shared']['logZ'],
        'dlnZ_free_anchored': out['M1_free_zT']['anchored']['logZ'] - out['M0a_LCDM']['anchored']['logZ'],
        'BF_zT5_shared':    float(np.exp(out['M1_zT5']['shared']['logZ']   - out['M0a_LCDM']['shared']['logZ'])),
        'BF_zT5_anchored':  float(np.exp(out['M1_zT5']['anchored']['logZ'] - out['M0a_LCDM']['anchored']['logZ'])),
        'BF_free_shared':   float(np.exp(out['M1_free_zT']['shared']['logZ']  - out['M0a_LCDM']['shared']['logZ'])),
        'BF_free_anchored': float(np.exp(out['M1_free_zT']['anchored']['logZ'] - out['M0a_LCDM']['anchored']['logZ'])),
    }
except Exception as e:
    out['nested_error'] = str(e)

out['interpretation'] = ('Controlled comparison using independent CEPH_DIST anchors and supplied covariance. Corpus sample/redshift convention retained; not an exact official-likelihood reproduction.')

Path('results').mkdir(exist_ok=True)
with open('results/audit_anchored_likelihood_corrected.json', 'w') as f:
    json.dump(out, f, indent=2)
print('Saved results/audit_anchored_likelihood_corrected.json')
