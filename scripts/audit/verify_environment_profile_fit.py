"""Independent bounded-intercept profiling and analytic-gradient G8 check."""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from scipy.stats import rankdata
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'steps'))
from step_03_11_env_modulated_kernel import load_pantheon,EnvKernelFit
from step_03_12_env_composite_kernel import local_density_rank
from verify_event_holdout import rankc
z,mb,mass,cal,cov,_=load_pantheon();fit=EnvKernelFit(z,mb,mass,cal,cov)
df=pd.read_csv('data/raw/pantheon_plus_shoes.dat',sep=r'\s+',comment='#')
mask=np.isfinite(df.zHD)&np.isfinite(df.m_b_corr)&(df.zHD>.001)&(df.zHD<3);df=df.loc[mask]
ambient=.5*(rankdata(mass)/len(z)+local_density_rank(z,df.RA.values,df.DEC.values,event_ids=df.CID.values));ambient-=ambient.mean()
ang=df.HOST_ANGSEP.values.astype(float);survey=df.IDSURVEY.values.astype(int)
ok=np.isfinite(ang)&(ang>0)&np.isin(survey,[1,10]);depth=np.zeros(len(z))
da=2*299792.458/70*(1-1/np.sqrt(1+z))/(1+z);raw=ang*da/(5*10**(.25*(mass-11)))
for s in np.unique(survey[ok]):
 ii=ok&(survey==s);r=rankdata(-np.log(np.maximum(raw[ii],1e-3)))/ii.sum();depth[ii]=r-r.mean()
quality=np.stack([rankc(df[c].values.astype(float)) for c in ['FITPROB','m_b_corr_err_DIAG','FITCHI2']])
ones=solve_triangular(fit.L,np.ones(len(z)),lower=True,check_finite=False);W=float(np.sum(ones**2))
records={}
for name,coords in [('quality',quality),('environment_quality',np.vstack([ambient,depth,quality]))]:
 def evaluate(theta,return_intercept=False):
  eps=theta[0];w=np.exp(theta[1:]@coords)
  c=fit.chi0+eps*w*fit.chi1
  mu=5*np.log10(c*(1+z))+25
  white=solve_triangular(fit.L,mb-mu,lower=True,check_finite=False)
  intercept=np.clip(np.sum(ones*white)/W,-21,-16)
  white-=intercept*ones
  precision=solve_triangular(fit.L.T,white,lower=False,check_finite=False)
  response=5/np.log(10)*w*fit.chi1/c
  jac=np.vstack([response,eps*response*coords])
  gradient=-2*np.sum(jac*precision,axis=1)
  value=float(np.sum(white**2))
  return intercept if return_intercept else (value,gradient)
 start=np.r_[.87,np.zeros(len(coords))]
 value,gradient=evaluate(start);h=1e-5
 finite=np.array([(evaluate(start+np.eye(len(start))[j]*h)[0]-evaluate(start-np.eye(len(start))[j]*h)[0])/(2*h) for j in range(len(start))])
 gradient_error=float(np.max(np.abs(gradient-finite))/max(1.,np.max(np.abs(gradient))));assert gradient_error<1e-6
 starts=[start];rng=np.random.default_rng(42)
 for _ in range(15):
  # Preserve the original multistart shape draws, consuming its M draw.
  starts.append(np.r_[rng.uniform(0,2),rng.uniform(-4,4,len(coords))]);rng.uniform(-19.8,-18.8)
 solutions=[minimize(evaluate,s,jac=True,method='L-BFGS-B',bounds=[(0,2)]+[(-4,4)]*len(coords),options={'ftol':1e-13,'gtol':1e-7,'maxiter':2000}) for s in starts]
 best=min(solutions,key=lambda r:r.fun);M=float(evaluate(best.x,True))
 w=np.exp(best.x[1:]@coords);direct=fit.chi2_given_w(best.x[0],M,w)
 assert abs(direct-best.fun)<1e-8
 records[name]=dict(chi2=float(best.fun),parameters=list(map(float,best.x))+[M],success=bool(best.success),message=str(best.message),analytic_gradient_relative_error=gradient_error,direct_likelihood_error=float(direct-best.fun),n_starts=len(starts),n_successful_starts=sum(r.success for r in solutions))
 print(name,records[name],flush=True)
records['environment_gain']=records['quality']['chi2']-records['environment_quality']['chi2']
Path('results/environment_profile_fit.json').write_text(json.dumps(records,indent=2)+'\n')
print('Gain',records['environment_gain'],flush=True)
