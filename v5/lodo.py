import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize, brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open('hse_bias.py').read().split('print("Allen+2008')[0])       # groups (grows,gresid), Allen arrays, Mstar_tot, fsolve, Hz, constants
    cd={}; exec(open('kernel_scenarios.py').read().split("# ---------- check NFW R500")[0],cd)
from mantz_data import MANTZ
PHI=.7; XC0=.4
# ---- X-COP ----
Dx=cd['build'](radii_mpc=[.5,1,1.5]); GX=cd['G']; xX=Dx['fb']*Dx['g']/Dx['gd']; densX=(Dx['g']/(GX*Dx['M5']/Dx['R5']**2))/Dx['w']; xiX=Dx['xi']
def xcop_res(p,idx=slice(None)):
    s=np.clip(densX,1e-6,None)**p; xs=xX/s; K=1+lam*np.exp(-xs)/(-np.expm1(-xs))**(.5+A*xiX); return np.log10(K*Dx['fb'])[idx]
# ---- groups ----
xiG=np.array([r['xi'] for r in grows]); stG=np.array([PHI*Mstar_tot(r['M5']/Msun)*Msun/(r['f25']*r['M25']) for r in grows])
def grp_res(p,idx=slice(None)): return np.array([gresid(r,stG[i],1,p=p[i]) for i,r in enumerate(grows)])[idx]
# ---- Allen ----
stA=PHI*Mstar_tot(M5/Msun)*Msun/(fo*M25); fbA=fo*(1+stA); xiA=xi.copy(); zA=z.copy(); lMA=lM.copy(); gdA=np.array([c_*Hz(zz)/(2*np.pi) for zz in zA])
def all_res(p,idx=slice(None)):
    pr=np.array([fsolve(m,r,x,g_,pp) for m,r,x,g_,pp in zip(M25,r25,xiA,gdA,p)]); return np.log(fbA/pr)[idx]
# ---- Mantz (held-out only, shape) ----
mz=[];
for n,zz,R,M,f,e in MANTZ:
    r_=R*kp_ if False else R*3.0857e19; m_=M*1e14*Msun; rc0=3*Hz(zz)**2/(8*np.pi*G_)
    fR=lambda l: 500*rc0*4*np.pi/3*(10**l)**3*mf(3.0*r_/10**l)/mf(3.0)-m_
    R5=10**brentq(fR,np.log10(r_*1.2),np.log10(r_*6)); M5m=500*rc0*4*np.pi/3*R5**3
    mz.append((zz,r_,m_,G_*M5m/R5**2/(c_*Hz(zz)/(2*np.pi)),c_*Hz(zz)/(2*np.pi),f,PHI*Mstar_tot(M5m/Msun)*Msun/(f*m_)))
mz=np.array(mz); zM,rM,mM,xiM,gdM,fM,stM=mz.T; lMM=np.log10(mM/Msun)
def man_res(p):
    pr=np.array([fsolve(m,r,x,g_,pp) for m,r,x,g_,pp in zip(mM,rM,xiM,gdM,p)]); return np.log(fM*(1+stM)/pr)
def P(model,th,xi_):
    return np.full_like(np.asarray(xi_,float),th[0]) if model=="const" else th[0]*(np.asarray(xi_,float)/XC0)**th[1]
def nll(th,model,use,bx=None,bg=None,ba=None):
    if th[0]<=0.03: return 1e9
    t=0
    if "X" in use: r=xcop_res(P(model,th,xiX),bx); t+=r.size*np.log(np.mean(r**2))
    if "G" in use: r=grp_res(P(model,th,xiG),bg); t+=r.size*np.log(np.mean(r**2))
    if "A" in use: r=all_res(P(model,th,xiA),ba); t+=r.size*np.log(np.mean(r**2))
    return t if np.isfinite(t) else 1e9
def fit(model,use,**kw):
    best=None
    for s in (([.5],) if model=="const" else ([.5,0.],[.5,.5],[.5,-.3])):
        r=minimize(nll,list(s),args=(model,use),method="Nelder-Mead",options=dict(xatol=2e-3,fatol=2e-3,maxiter=250))
        if best is None or r.fun<best.fun: best=r
    return best
def slope(y,x): return np.polyfit(x,y,1)[0]
def evaluate(model,th):
    out={}
    rX=xcop_res(P(model,th,xiX)); out['X-COP rms']=np.sqrt(np.mean(rX**2))
    rG=grp_res(P(model,th,xiG)); out['groups rms/mean']="%.3f/%+.3f"%(np.sqrt(np.mean(rG**2)),rG.mean())
    rA=all_res(P(model,th,xiA)); pr=fbA/np.exp(rA); out['Allen rms(ln)/mean']="%.3f/%+.3f"%(np.sqrt(np.mean(rA**2)),rA.mean()); out['Allen mass slope model']=slope(np.log(pr),lMA)
    rM_=man_res(P(model,th,xiM)); d=rM_-rM_.mean(); fp=fM*(1+stM)/np.exp(rM_); out['Mantz shape rms']=np.sqrt(np.mean(d**2)); out['Mantz mass slope model']=slope(np.log(fp),lMM)
    return out
print("observed mass slopes: Allen %+.2f ; Mantz %+.2f  (stars phi=0.7, Kravtsov relation)"%(slope(np.log(fbA),lMA),slope(np.log(fM*(1+stM)),lMM)))
SC=[("train groups+X-COP -> predict Allen,Mantz",("G","X")),("train Allen+X-COP -> predict groups,Mantz",("A","X")),("train groups+Allen -> predict X-COP,Mantz",("G","A")),("train all three -> Mantz (held out)",("G","A","X"))]
for nm,use in SC:
    print("\n##",nm)
    for model in ("const","pow"):
        r=fit(model,use); ev=evaluate(model,r.x)
        held=[k for k in ev if (("Allen" in k and "A" not in use) or ("groups" in k and "G" not in use) or ("X-COP" in k and "X" not in use) or "Mantz" in k)]
        print("  %-5s params %s | -2lnL(train) %.1f | held-out: %s"%(model,np.round(r.x,2),r.fun,"; ".join("%s=%s"%(k,("%.3f"%ev[k]) if isinstance(ev[k],float) else ev[k]) for k in held)))
print("\n## bootstrap of (p0,q) with all three datasets (resampling systems within each dataset), 120 reps")
rng=np.random.default_rng(5); qs=[];p0s=[]
for _ in range(120):
    bx=rng.integers(0,xiX.shape[0],xiX.shape[0]); bg=rng.integers(0,len(grows),len(grows)); ba=rng.integers(0,len(zA),len(zA))
    f=lambda th: nll(th,"pow",("X","G","A"),bx,bg,ba)
    r=minimize(f,[.47,.36],method="Nelder-Mead",options=dict(xatol=3e-3,fatol=3e-3,maxiter=120)); p0s.append(r.x[0]); qs.append(r.x[1])
qs=np.array(qs); print(" q: median %.2f, 68%% [%.2f, %.2f], 95%% [%.2f, %.2f]; P(q>0)=%.3f ; p0 median %.2f [%.2f, %.2f]"%(np.median(qs),*np.percentile(qs,[16,84]),*np.percentile(qs,[2.5,97.5]),np.mean(qs>0),np.median(p0s),*np.percentile(p0s,[16,84])))
