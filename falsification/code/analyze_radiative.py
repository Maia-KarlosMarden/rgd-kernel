import os as _os; _os.chdir(_os.path.dirname(_os.path.abspath(__file__)))  # paths relativos a v5/code
import numpy as np, math, json
from scipy.optimize import minimize_scalar
from sklearn.model_selection import KFold
from sklearn.linear_model import Ridge
G=6.67430e-11; KPC=3.085677581491367e19; c=299792458.
meta={}
for line in open('data/SPARC_Lelli2016c.mrt'):
 if not line.startswith('   '):continue
 try:
  name=line[0:11].strip(); vals=line[11:].split()
  if len(vals)<16:continue
  T,D,eD,fd,inc,einc,L,eL,Re,SBef,Rd,SBd,MHI,RHI,Vflat,eVflat,Q=map(float,vals[:17])
  if name: meta[name]=dict(L=L,Re=Re,Rd=Rd,SBd=SBd,MHI=MHI,Vflat=Vflat,Q=Q)
 except: pass
print('metadata',len(meta),list(meta)[:5])
points=[]
for line in open('data/SPARC_Table2.txt'):
 if not line or not line[:11].strip() in meta:continue
 try:
  name=line[:11].strip(); vals=list(map(float,line[11:].split())); D,r,v,ev,vg,vd,vbul,sbd,sbb=vals
  m=meta[name]
  if m['Q']>2 or r<=0 or v<=0 or m['L']<=0 or m['Rd']<=0:continue
  vb2=vg*abs(vg)+.5*vd**2+.7*vbul**2
  if vb2<=0:continue
  gb=(vb2*1e6)/(r*KPC); go=(v*v*1e6)/(r*KPC)
  if gb<=0 or go<=0:continue
  x=r/m['Rd']; cfrac=(1-(1+x)*math.exp(-x)); cfrac=max(cfrac,1e-8)
  # 3 proxies: concentration from disk radial light enclosed; geometrical deflection proxy; 3.6um effective irradiance proxy
  C=m['SBd'] * cfrac/(x*x)  # projected average luminosity surface density proxy (constant factors irrelevant)
  geom=4*vb2*1e6/c**2
  flux=m['L']*1e9/(4*math.pi*r*r) # Lsun(3.6) per kpc^2
  points.append(dict(name=name,r=r,gb=gb,go=go,C=C,G=geom,L=flux,lum=m['L'],v=v,vb=math.sqrt(vb2),Q=m['Q']))
 except Exception as e: pass
names=sorted(set(p['name'] for p in points)); print('points',len(points),'galaxies',len(names))
# galaxies with >=5
names=[n for n in names if sum(p['name']==n for p in points)>=5]
points=[p for p in points if p['name'] in set(names)]
print('filtered',len(points),len(names))
y=np.log10([p['go'] for p in points]); gb=np.array([p['gb'] for p in points]); arr=np.array([[math.log10(p[k]) for k in ['C','G','L']] for p in points]); center=np.median(arr,axis=0); scale=np.std(arr,axis=0); Z=(arr-center)/scale
# triad factor (normalized product) = log C +log G+log L, and multi feature (3 linear + interactions)
X1=Z.sum(axis=1).reshape(-1,1)/math.sqrt(3)
X3=np.column_stack([Z, Z[:,0]*Z[:,1]*Z[:,2]])
nameidx=np.array([names.index(p['name']) for p in points]); folds=list(KFold(n_splits=5,shuffle=True,random_state=42).split(names))
def pred(a,g):return np.log10(g/(1-np.exp(-np.sqrt(g/a))))
models=['RAR','RAR+produto_CGL','RAR+componentes_interacao']
out={m:np.full(len(points),np.nan) for m in models}; fitted=[]
for fold,(itr,ite) in enumerate(folds):
 tr=np.isin(nameidx,itr); te=np.isin(nameidx,ite)
 def cost(loga):
  er=pred(10**loga,gb[tr])-y[tr]
  return np.mean([np.mean(er[nameidx[tr]==j]**2) for j in itr if np.any(nameidx[tr]==j)])
 opt=minimize_scalar(cost,bounds=(-12,-9),method='bounded',options={'xatol':1e-9}); a=10**opt.x
 base_tr=pred(a,gb[tr]); base_te=pred(a,gb[te]);out['RAR'][te]=base_te
 for model,X in [('RAR+produto_CGL',X1),('RAR+componentes_interacao',X3)]:
  # fixed regularization, galaxy equal weights
  counts=np.bincount(nameidx[tr],minlength=len(names)); weights=1/counts[nameidx[tr]]
  clf=Ridge(alpha=10.0,fit_intercept=True)
  clf.fit(X[tr],y[tr]-base_tr,sample_weight=weights)
  out[model][te]=base_te+clf.predict(X[te])
 fitted.append((fold,a, len(itr),sum(te)))
print('folds',fitted)
for m in models:
 er=out[m]-y
 per=[np.mean(er[nameidx==i]**2) for i in range(len(names))]
 print(m,'rmse dex',np.sqrt(np.mean(per)),'all point rmse',np.sqrt(np.mean(er**2)))
# full fit for scatter, BTFR summary independent
opt=minimize_scalar(lambda la:np.mean([np.mean((pred(10**la,gb[nameidx==i])-y[nameidx==i])**2) for i in range(len(names))]),bounds=(-12,-9),method='bounded'); afull=10**opt.x
btfr=[]
for n in names:
 m=meta[n]
 if m['Vflat']>0:
  mb=(0.5*m['L']+1.33*m['MHI'])*1e9*1.98847e30
  btfr.append((n,((m['Vflat']*1e3)**4/(G*mb))))
print('a full',afull,'BTFR',len(btfr),np.median([x[1] for x in btfr]))
result={'points':points,'meta':meta,'names':names,'folds':fitted,'metrics':{m:{'RMSE_per_galaxy':float(np.sqrt(np.mean([np.mean((out[m][nameidx==i]-y[nameidx==i])**2) for i in range(len(names))]))),'RMSE_points':float(np.sqrt(np.mean((out[m]-y)**2)))} for m in models},'a_full':afull,'BTFR_med':float(np.median([x[1] for x in btfr])),'BTFR':btfr,'sample':[]}
for n in ['NGC2403','NGC3198','NGC2841','NGC7331','NGC7814','DDO154']:
 if n in names:
  idx=np.where(nameidx==names.index(n))[0]; ix=idx[-1]
  result['sample'].append(dict(galaxy=n,R=points[ix]['r'],g_obs=points[ix]['go'],g_b=points[ix]['gb'],A=points[ix]['go']/points[ix]['gb']-1,ratio=points[ix]['go']/points[ix]['gb'],RAR=10**out['RAR'][ix],Triple=10**out['RAR+componentes_interacao'][ix]))
with open('data/SPARC_tripla_test_results.json','w') as f:json.dump(result,f,default=lambda x: x.item() if hasattr(x,'item') else str(x))
