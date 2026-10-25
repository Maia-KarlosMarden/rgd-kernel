"""PRE-REGISTERED (written before running). Is the radial scale factor s a function of a dimensionless ratio?
Analytic identity checked numerically: with Delta500 as reference, s = dd^(1/2) means g_eff = g_dagger*dd^(1/2) = c*Omega_dyn/(2*pi*sqrt(250)),
 Omega_dyn = sqrt(G M(<r)/r^3): H(z) cancels -> the threshold follows the LOCAL dynamical frequency instead of the Hubble rate (p=1/2 is the unique power giving a frequency).
Test (direct fit of ln f_bar, not of p_i; same 254+ systems as ref51): ln s = F(w) for dimensionless ratios w in
  dd=Delta/Delta500, gr=g(r)/g500, mr=M(<r)/M500, rr=r/R500, gg=g(r)/g_dagger, xi=g500/g_dagger.
Models: (a) s=1; (b) s=dd^0.5 (archived); (c) s=w^q free q for each w; (d) s=exp(c)dd^0.5 (c <-> free reference overdensity Delta_ref=500 exp(-2c));
 (e) s=dd^p exp(c); (f) s=dd^p xi^a ; (g) (c)+XMM offset in ln s.
Likelihood: ln f_obs ~ N(ln f_pred, sf_i^2 + tau^2).  Criteria: AIC (k params + tau) and leave-one-sample-out (LOSO) log score. Parameters never tuned on a held-out sample.
"""
import io,contextlib,sys,numpy as np,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import ref51 as Q
R=Q.R; D=R.D; Gc,cc,lam_,A_=R.Gc,R.cc,R.lam_,R.A_
from scipy.optimize import brentq, minimize
def gd_of(d): return cc*d['H0']*R.E(d['z'])/(2*np.pi)
def fpred_s(d,s):
    gd=gd_of(d); xi=Gc*d['M5']/d['R5']**2/gd; g=Gc*d['M']/d['r']**2; n=.5+A_*xi
    h=lambda f: f*(1+lam_*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**n)-1
    try: return brentq(h,1e-6,5.0)
    except ValueError: return np.nan
# ratios
for d in D:
    d['w']=dict(dd=d['dd'],gr=(d['M']/d['M5'])/(d['r']/d['R5'])**2,mr=d['M']/d['M5'],rr=d['r']/d['R5'],
                gg=Gc*d['M']/d['r']**2/gd_of(d),xi=Gc*d['M5']/d['R5']**2/gd_of(d))
    d['y']=np.log(d['fb'])
# identity check s=dd^0.5 <-> c*Omega/(2 pi sqrt250) (valid when M5,R5 obey Delta=500: test on XXL systems)
ch=[]
for d in D:
    if d['samp']!="XXL C1": continue
    Om=np.sqrt(Gc*d['M']/d['r']**3); lhs=gd_of(d)*np.sqrt(d['dd']); rhs=cc*Om/(2*np.pi*np.sqrt(250)); ch.append(lhs/rhs)
    if len(ch)>40: break
print("identity g_dag*sqrt(dd) / (c*Omega/(2pi sqrt250)) on 41 XXL systems: mean %.4f, spread %.4f (Omega^2=GM/r^3, uses E(z) with Om=0.3 for g_dag and Om=0.28 for R500)"%(np.mean(ch),np.std(ch)))
samp=np.array([d['samp'] for d in D]); xmm=np.array([d['instr'] for d in D],float); sf=np.array([d['sf'] for d in D]); y=np.array([d['y'] for d in D])
W={k:np.array([d['w'][k] for d in D]) for k in D[0]['w']}
def lns_fn(model,th,idx):
    if model=='none': return np.zeros(len(idx))
    if model=='dd05': return 0.5*np.log(W['dd'][idx])
    if model[:2]=='q_': return th[0]*np.log(W[model[2:]][idx])
    if model=='c05': return th[0]+0.5*np.log(W['dd'][idx])
    if model=='pc': return th[0]*np.log(W['dd'][idx])+th[1]
    if model=='pxi': return th[0]*np.log(W['dd'][idx])+th[1]*np.log(W['xi'][idx])
    if model=='c05x': return th[0]+0.5*np.log(W['dd'][idx])+th[1]*xmm[idx]
    if model=='px': return th[0]*np.log(W['dd'][idx])+th[1]*xmm[idx]
npar={'none':0,'dd05':0,'c05':1,'pc':2,'pxi':2,'c05x':2,'px':2}
def resid(model,th,idx):
    ls=lns_fn(model,th,idx); out=np.empty(len(idx))
    for j,i in enumerate(idx):
        fp=fpred_s(D[i],np.exp(ls[j])); out[j]=y[i]-np.log(fp) if np.isfinite(fp) and fp>0 else np.nan
    return out
def nll(model,th,lt,idx):
    r=resid(model,th,idx)
    if np.any(~np.isfinite(r)): return 1e9
    v=sf[idx]**2+np.exp(2*lt); return 0.5*np.sum(np.log(2*np.pi*v)+r**2/v)
def fit(model,idx,th0):
    k=npar.get(model,1)
    f=lambda x: nll(model,x[:k],x[k],idx)
    o=minimize(f,np.r_[th0[:k],np.log(0.15)],method='Nelder-Mead',options=dict(xatol=1e-3,fatol=1e-3,maxiter=400))
    return o.x[:k],np.exp(o.x[k]),o.fun
allidx=np.arange(len(D)); samples=sorted(set(samp))
if __name__=='__main__':
    models=[('none',[]),('dd05',[]),('c05',[0.]),('pc',[.5,0.]),('pxi',[.5,0.]),('px',[.5,0.]),('c05x',[0.,0.])]+[('q_'+k,[.3]) for k in ('dd','gr','mr','rr','gg','xi')]
    print("N systems:",len(D))
    print("\n%-8s k  tau    nll     AIC    dAIC  LOSO(sum)  dLOSO  params"%"model")
    rows=[]
    for m,th0 in models:
        th,tau,v=fit(m,allidx,th0); k=npar.get(m,1); aic=2*v+2*(k+1)
        loso=0
        for sm in samples:
            te=np.where(samp==sm)[0]; tr=np.where(samp!=sm)[0]
            tht,taut,_=fit(m,tr,th0)
            loso-=nll(m,tht,np.log(taut),te)
        rows.append((m,k,tau,v,aic,loso,th))
    a0=rows[0][4]; l0=rows[0][5]
    for m,k,tau,v,aic,loso,th in rows:
        print("%-8s %d  %.3f  %7.1f %7.1f %+7.1f  %8.1f  %+7.1f  %s"%(m,k,tau,v,aic,aic-a0,loso,loso-l0,np.round(th,3)))
