import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
import numpy as np, json, math
from scipy.optimize import minimize, minimize_scalar
from sklearn.model_selection import KFold
from collections import Counter
D=json.load(open('data/SPARC_tripla_test_results.json'))
pts=D['points']; names=D['names']; n=len(pts)
idx=np.array([names.index(p['name']) for p in pts]); gb=np.array([p['gb'] for p in pts]); go=np.array([p['go'] for p in pts]); y=np.log10(go)
X=np.log10(np.array([[p['C'],p['G'],p['L']] for p in pts])); cts=np.bincount(idx,minlength=len(names)); wt=1/cts[idx]; wt/=wt.sum()
folds=list(KFold(n_splits=5,shuffle=True,random_state=42).split(names))
def rar(a, gg): return gg/(-np.expm1(-np.sqrt(gg/a)))
def fit_a(mask):
  def obj(z):
    e=np.log10(rar(10**z,gb[mask]))-y[mask]; ww=wt[mask];return np.sum(ww*e*e)/ww.sum()
  return 10**minimize_scalar(obj,bounds=(-12,-9),method='bounded',options={'xatol':1e-8}).x

def fit(theta0, xx, mask, gbbase, pat, ridge):
    def calc(t):
      z=np.clip(t[0]+xx[mask]@t[1:],-6,6)
      pred=gb[mask]+gbbase[mask]*np.power(10,z)
      err=np.log10(pred)-y[mask]
      val=np.sum(wt[mask]*err**2)/wt[mask].sum()
      return val+ridge*np.sum(t[1:]**2)
    bounds=[(-2,2)]+[(0,3) if a=='+' else (-3,0) for a in pat]
    res=minimize(calc,theta0,method='L-BFGS-B',bounds=bounds,options={'maxiter':180,'ftol':1e-11})
    return res.x
patterns=[''.join(v) for v in __import__('itertools').product('+-', repeat=3)]
preds={k:np.full(n,np.nan) for k in ['RAR','free']+patterns}; coeffs={k:[] for k in ['free']+patterns}; a_folds=[]
for k,(itr,ite) in enumerate(folds):
 tr=np.isin(idx,itr); te=np.isin(idx,ite)
 a=fit_a(tr);a_folds.append(a)
 bbase=rar(a,gb)-gb
 center=np.median(X[tr],axis=0);xx=X-center
 preds['RAR'][te]=rar(a,gb[te])
 for pat in ['free']+patterns:
  if pat=='free':
   bounds=[(-2,2)]+[(-3,3)]*3
   def calc(t):
      pp=gb[tr]+bbase[tr]*10**np.clip(t[0]+xx[tr]@t[1:],-6,6)
      er=np.log10(pp)-y[tr]
      return np.sum(wt[tr]*er**2)/wt[tr].sum()+0.002*np.sum(t[1:]**2)
   res=minimize(calc,np.zeros(4),method='L-BFGS-B',bounds=bounds,options={'maxiter':180,'ftol':1e-11});th=res.x
  else: th=fit(np.zeros(4),xx,tr,bbase,pat,0.002)
  pp=gb[te]+bbase[te]*10**np.clip(th[0]+xx[te]@th[1:],-6,6)
  preds[pat][te]=pp
  coeffs[pat].append({'fold':k+1,'a0':a,'log10K':float(th[0]),'p':float(th[1]),'q':float(th[2]),'s':float(th[3])})
 print('completed fold',k+1,flush=True)

def metrics(pp):
 er=np.log10(pp)-y
 return {'RMSE_galaxy':float(np.sqrt(np.mean([np.mean(er[idx==i]**2) for i in range(len(names))]))),
         'RMSE_points':float(np.sqrt(np.mean(er*er))),'pergalaxy':list(map(float,[np.sqrt(np.mean(er[idx==i]**2)) for i in range(len(names))]))}
mets={k:metrics(p) for k,p in preds.items()}
for k in sorted(mets,key=lambda k:mets[k]['RMSE_galaxy']): print(k,round(mets[k]['RMSE_galaxy'],6), 'coeff',np.round(np.mean([[q['p'],q['q'],q['s']] for q in coeffs.get(k,[])],axis=0),4) if k!='RAR' else '',flush=True)
print('corr',np.corrcoef(X.T));print('fraction negative',np.mean(go<=gb),np.sum(go<=gb))
result={'metrics':mets,'coefficients':coeffs,'a_folds':a_folds,'corr':np.corrcoef(X.T).tolist(),'data': {'num_galaxies':len(names),'n_points':n,'neg_residual':int(np.sum(go<=gb))},'method':'Fold grouped by galaxy, a0 refit in training, multiplicative correction to positive RAR residual; L-BFGS-B sign constrained; regularization 0.002*sum(exponent^2); train-centered proxy log10; p,q,s bounds ±3; gobs log10 loss per galaxy; 5 folds random_state 42', 'galaxies':names}
json.dump(result,open('data/SPARC_expoentes_livres_resultados.json','w'))
