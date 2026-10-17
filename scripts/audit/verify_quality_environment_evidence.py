"""Independent quality-controlled evidence under the existing +-1 gate priors."""
import runpy,json
from pathlib import Path
import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import ndtr,logsumexp
from scipy.stats import qmc,norm,chi2,multivariate_t
state=runpy.run_path('scripts/audit/verify_environment_profile_fit.py')
fit=state['fit'];mb=state['mb'];z=state['z'];ones=state['ones'];W=state['W']
quality=state['quality'];ambient=state['ambient'];depth=state['depth']
modes=state['records'];sigma=1/np.sqrt(W)
checks=[]
for name,coordinates in [('quality',quality),('environment_quality',np.vstack([ambient,depth,quality]))]:
 center=np.array(modes[name]['parameters'][:-1]);dim=len(center)
 def marginal(theta):
  theta=np.atleast_2d(theta);eps=theta[:,0]
  exponents=np.sum(coordinates[:,:,None]*theta[:,1:].T[:,None,:],axis=0)
  weight=np.exp(exponents)
  distance=(fit.chi0[:,None]+fit.chi1[:,None]*eps*weight)*(1+z[:,None])
  mu=5*np.log10(distance)+25
  white=solve_triangular(fit.L,mb[:,None]-mu,lower=True,check_finite=False)
  mhat=np.sum(ones[:,None]*white,axis=0)/W
  chi=np.sum((white-ones[:,None]*mhat)**2,axis=0)
  prob=ndtr((-16-mhat)/sigma)-ndtr((-21-mhat)/sigma)
  return -.5*chi+.5*np.log(2*np.pi/W)-np.log(5)+np.log(prob)
 f=lambda x:-float(marginal(x)[0]);f0=f(center);h=1e-4
 H=np.zeros((dim,dim))
 for i in range(dim):
  ei=np.eye(dim)[i]*h;H[i,i]=(f(center+ei)-2*f0+f(center-ei))/h**2
  for j in range(i):
   ej=np.eye(dim)[j]*h
   H[i,j]=H[j,i]=(f(center+ei+ej)-f(center+ei-ej)-f(center-ei+ej)+f(center-ei-ej))/(4*h*h)
 assert np.linalg.eigvalsh(H).min()>0
 shape=np.linalg.inv(H);chol=np.linalg.cholesky(shape)
 for seed in (101,202):
  u=qmc.Sobol(dim+1,scramble=True,seed=seed).random_base2(14)
  gaussian=norm.ppf(u[:,:dim])
  shifted=np.sum(gaussian[:,None,:]*chol[None,:,:],axis=2)
  sample=center+shifted/np.sqrt(chi2.ppf(u[:,dim],5)/5)[:,None]
  assert np.all(np.isfinite(sample))
  proposal=multivariate_t.logpdf(sample,loc=center,shape=shape,df=5)
  prior_ok=(sample[:,0]>=0)&(sample[:,0]<=2)&np.all(np.abs(sample[:,1:])<=1,axis=1)
  ll=np.full(len(sample),-np.inf)
  for start in range(0,len(sample),128):
   sl=slice(start,start+128);active=prior_ok[sl]
   if active.any():
    sub=ll[sl];sub[active]=marginal(sample[sl][active]);ll[sl]=sub
  lw=ll-proposal
  logz=logsumexp(lw)-np.log(len(sample))-dim*np.log(2)
  ess=np.exp(2*logsumexp(lw)-logsumexp(2*lw))
  record=dict(model=name,seed=seed,n=len(sample),logZ=float(logz),importance_ESS=float(ess),prior_accepted_samples=int(prior_ok.sum()),posterior_scales=np.sqrt(np.diag(shape)).tolist())
  checks.append(record);print(record,flush=True)
ratios=[]
for seed in (101,202):
 values={x['model']:x['logZ'] for x in checks if x['seed']==seed}
 ratios.append(dict(seed=seed,BF_environment_extension_vs_quality=float(np.exp(values['environment_quality']-values['quality']))))
out=dict(scope='Same G8 likelihood and priors; independent bounded-intercept integral and importance estimator',priors=dict(epsilon=[0,2],all_gate_coefficients=[-1,1],magnitude=[-21,-16]),checks=checks,comparisons=ratios)
Path('results/quality_environment_evidence_independent.json').write_text(json.dumps(out,indent=2)+'\n');print(ratios,flush=True)
