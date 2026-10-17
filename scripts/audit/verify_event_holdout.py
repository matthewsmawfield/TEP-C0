"""Matched random holdout with physical events kept in a single fold.

Coordinates retain the specified full-catalog, outcome-free construction.
This checks parameter transfer within this catalog, not external replication.
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'steps'))
from step_03_12_env_composite_kernel import local_density_rank, EnvKernelFit, load_pantheon

def rankc(values):
    ok = np.isfinite(values)
    result = np.zeros(len(values))
    ranks = rankdata(values[ok])/ok.sum()
    result[ok] = ranks-ranks.mean()
    return result

def main():
    df = pd.read_csv('data/raw/pantheon_plus_shoes.dat', sep=r'\s+', comment='#')
    valid = np.isfinite(df.zHD)&np.isfinite(df.m_b_corr)&(df.zHD>.001)&(df.zHD<3.)
    df = df[valid]
    z, mb, mass, cal, cov, _ = load_pantheon()
    ids = df.CID.to_numpy()
    density = local_density_rank(z, df.RA.to_numpy(), df.DEC.to_numpy(), event_ids=ids)
    ambient = .5*(rankdata(mass)/len(mass)+density)
    ambient -= ambient.mean()
    ang = df.HOST_ANGSEP.to_numpy(float)
    survey = df.IDSURVEY.to_numpy(int)
    ok = np.isfinite(ang)&(ang>0)&np.isin(survey,[1,10])
    da = (2*299792.458/70)*(1-1/np.sqrt(1+z))/(1+z)
    raw = ang*da/(5*10**(.25*(mass-11)))
    depth = np.zeros(len(z))
    for s in np.unique(survey[ok]):
        mask = ok&(survey==s)
        values = rankdata(-np.log(np.maximum(raw[mask],1e-3)))/mask.sum()
        depth[mask] = values-values.mean()
    quality = np.stack([rankc(df[c].to_numpy(float)) for c in ['FITPROB','m_b_corr_err_DIAG','FITCHI2']])
    coordinates = np.vstack([ambient,depth,quality])
    events, inverse = np.unique(ids, return_inverse=True)
    event_fold = np.zeros(len(events),int)
    permutation = np.random.default_rng(3).permutation(len(events))
    for k in range(5):
        event_fold[permutation[k::5]] = k
    folds = event_fold[inverse]
    records = []
    for k in range(5):
        test = folds==k
        train = ~test
        assert not (set(ids[test])&set(ids[train]))
        fit = EnvKernelFit(z[train],mb[train],mass[train],cal[train],cov[np.ix_(train,train)])
        models = []
        test_fit = EnvKernelFit(z[test],mb[test],mass[test],cal[test],cov[np.ix_(test,test)])
        for coords in (quality,coordinates):
            result = fit.fit(lambda p: np.exp(p@coords[:,train]), [(-4.,4.)]*len(coords))
            w = np.exp(result.x[1:-1]@coords[:,test])
            residual = mb[test]-test_fit.mu_model(result.x[0],w)-result.x[-1]
            white = solve_triangular(test_fit.L,residual,lower=True)
            models.append(dict(training_chi2=float(result.fun),parameters=result.x.tolist(),held_out_chi2=float(np.sum(white*white))))
        row = dict(fold=k,held_out_rows=int(test.sum()),held_out_events=int(sum(event_fold==k)),shared_events=0,
                   quality_only=models[0],environment_quality=models[1],held_out_gain=models[0]['held_out_chi2']-models[1]['held_out_chi2'])
        records.append(row)
        print(json.dumps(row),flush=True)
    output = dict(seed=3,partition_unit='CID physical event',rows=len(z),events=len(events),
                  coordinate_scope='Outcome-free full catalog coordinates; unchanged construction',
                  evaluation='Frozen fitted parameters; test marginal covariance; no test outcome refit',
                  folds=records,total_held_out_gain=sum(r['held_out_gain'] for r in records))
    Path('results/event_grouped_holdout.json').write_text(json.dumps(output,indent=2)+'\n')
    print('TOTAL',output['total_held_out_gain'],flush=True)

if __name__ == '__main__':
    main()
