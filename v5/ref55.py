"""PRE-REGISTERED (written before running).  Halo-level law F(xi): baryon fraction inside r500 vs xi = g500/g_dagger, s=1 at r500 so radial law is irrelevant.
f_b solves f*K(f*xi; n(xi)) = 1,  K = 1 + lam exp(-x)/(1-exp(-x))^n.
Data (stars: Kravtsov, phi=1): HSC-XXL population relation at z=0.3 (M=2e13,3e13,5.6e13,1e14,2e14,3e14; MC errors, treated as independent), Vikhlinin06 f500 (10), Ettori19 f500 (9),
 groups f500 (Lovisari15 XMM + Sun09 Chandra, 43), XXL-365 C1 per-system IF f(<500kpc)=f(<r500) (flagged, separate), eFEDS bins (flagged, separate).
Families: N0 archived (lam=.999, n=.5+.4068 xi); N1 refit A; N2 n=const refit; N3 refit n0,A; N4 refit lam,n0,A.  Hierarchical scatter tau on ln f_b, per-sample offset NOT allowed in the main fit.
Criteria: AIC and leave-one-sample-out (LOSO) log score.  Main fit sample set = {HSC-XXL, Vikhlinin, Ettori, Groups-L, Groups-S}; XXL-C1 and eFEDS reported as external."""
import io,contextlib,sys,os,numpy as np,warnings; warnings.filterwarnings("ignore")
STARS=os.environ.get('STARS','krav')   # declared sensitivity: 'xxl' = HSC-XXL measured stellar fraction instead of Kravtsov
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    exec(open('stars_test.py').read())
    GR_=list(grows)
    import blind_xxl_efeds as B
    import predict_frozen as pf
from scipy.optimize import brentq, minimize
Mpc_=3.0857e22; Ms_=1.989e30
E=lambda z: np.sqrt(.3*(1+z)**3+.7)
def stfrac(M5,z):
    return Mstar_tot(M5)/M5 if STARS=='krav' else B.xxl_frac(M5,z)[1]
rows=[]   # (sample, ln f_obs, sigma_ln, M500[Msun], z, xi)
def xi_of(M5,z):   # M5 in Msun
    H=70e3/Mpc_*E(z); rho=3*H**2/(8*np.pi*6.674e-11); R5=(3*M5*Ms_/(4*np.pi*500*rho))**(1/3)
    return 6.674e-11*M5*Ms_/R5**2/(3e8*H/(2*np.pi))
rng=np.random.default_rng(5)
for M in (2e13,3e13,5.6e13,1e14,2e14,3e14):
    g,s=B.xxl_frac(M,.3,rng); fb=g+s; rows.append(("HSC-XXL",np.log(np.median(fb)),np.std(np.log(fb)),M,.3,xi_of(M,.3)))
Hv=lambda z: 71e3/Mpc_*E(z)
V=[("A133",.0569,998,1.24,3.14,.065,.084),("A383",.1883,956,1.68,3.10,.090,.122),("A478",.0881,1359,4.23,7.83,.096,.118),("A907",.1603,1117,2.30,4.71,.090,.121),("A1413",.1429,1339,3.08,7.78,.092,.104),("A1795",.0622,1283,2.66,6.57,.088,.098),("A1991",.0592,753,.60,1.28,.068,.099),("A2029",.0779,1380,4.41,8.29,.090,.121),("A2390",.2302,1448,3.50,10.88,.127,.138),("MKW4",.0199,635,.30,.74,.044,.063)]
for n,z,r5,M25,M5,f25,f5 in V:
    st=stfrac(M5*1e14,z); rows.append(("Vikhlinin",np.log(f5+st),0.10,M5*1e14,z,xi_of(M5*1e14,z)))
Xe=[("A1644",.0473,3.48,6.69,.128,.126),("A1795",.0622,4.63,6.53,.139,.144),("A2029",.0773,8.65,12.25,.141,.152),("A2142",.0909,8.95,13.64,.158,.168),("A2255",.0809,5.26,10.33,.153,.146),("A3158",.0597,4.26,6.63,.145,.155),("A644",.0704,5.66,7.67,.132,.139),("A85",.0555,5.65,8.50,.150,.159),("RXC1825",.0650,4.08,6.15,.133,.155)]
for n,z,M5,M2,f5,f2 in Xe:
    st=stfrac(M5*1e14,z); rows.append(("Ettori",np.log(f5+st),0.10,M5*1e14,z,xi_of(M5*1e14,z)))
for r in GR_:
    M5=r['M5']/Msun; st=stfrac(M5,r['z']); rows.append(("Groups-L" if r['src']=="L" else "Groups-S",np.log(r['f5']+st),0.10,M5,r['z'],xi_of(M5,r['z'])))
samp=np.array([r[0] for r in rows]); y=np.array([r[1] for r in rows]); sg=np.array([r[2] for r in rows]); M=np.array([r[3] for r in rows]); z=np.array([r[4] for r in rows]); xi=np.array([r[5] for r in rows])
print("N=%d"%len(y)," ".join("%s:%d"%(s,(samp==s).sum()) for s in sorted(set(samp))),"  xi range %.3f..%.3f"%(xi.min(),xi.max()))
def fF(xi_,lam,n0,A):
    out=np.empty(len(xi_))
    for j,x in enumerate(xi_):
        n=n0+A*x; h=lambda f: f*(1+lam*np.exp(-f*x)/(-np.expm1(-f*x))**n)-1
        try: out[j]=brentq(h,1e-6,5.0)
        except ValueError: out[j]=np.nan
    return out
def nll(th,mode,idx,xi_,y_,s_):
    lam,n0,A=th[:3]
    f=fF(xi_[idx],lam,n0,A)
    if np.any(~np.isfinite(f)): return 1e9
    r=y_[idx]-np.log(f); v=s_[idx]**2+np.exp(2*th[3])
    return 0.5*np.sum(np.log(2*np.pi*v)+r**2/v)
def full(th,mode):
    lam=.999; n0=.5; A=.4068
    if mode=='N1': A=th[0]
    elif mode=='N2': n0=th[0]; A=0
    elif mode=='N3': n0,A=th[0],th[1]
    elif mode=='N4': lam,n0,A=th[0],th[1],th[2]
    return lam,n0,A
spec={'N0':(0,[]),'N1':(1,[.4]),'N2':(1,[.6]),'N3':(2,[.5,.4]),'N4':(3,[.999,.5,.4])}
def fitm(mode,idx):
    k,t0=spec[mode]
    if k==0: 
        f=lambda lt: nll(list(full([],mode))+[lt[0]],mode,idx,xi,y,sg); o=minimize(f,[np.log(.1)],method='Nelder-Mead'); return [],np.exp(o.x[0]),o.fun
    def f(x):
        th=full(x[:k],mode)
        if mode=='N4' and not (0.05<x[0]<5): return 1e9
        return nll(list(th)+[x[k]],mode,idx,xi,y,sg)
    o=minimize(f,np.r_[t0,np.log(.1)],method='Nelder-Mead',options=dict(maxiter=800,xatol=1e-4,fatol=1e-4)); return list(o.x[:k]),np.exp(o.x[k]),o.fun
allidx=np.arange(len(y)); samples=sorted(set(samp))
print("\nmodel  k  tau    nll     AIC    LOSO   params")
res0=None
for mode in ('N0','N1','N2','N3','N4'):
    th,tau,v=fitm(mode,allidx); k=spec[mode][0]; aic=2*v+2*(k+1); lo=0
    for sm in samples:
        te=np.where(samp==sm)[0]; tr=np.where(samp!=sm)[0]; tht,taut,_=fitm(mode,tr)
        lo-=nll(list(full(tht,mode))+[np.log(taut)],mode,te,xi,y,sg)
    print("%-5s %d  %.3f  %6.1f  %6.1f  %6.1f  %s"%(mode,k,tau,v,aic,lo,np.round(th,3)))
    if mode=='N0': res0=(th,tau)
# per-sample residuals under N0 (median ratio obs/pred)
lam,n0,A=.999,.5,.4068; f=fF(xi,lam,n0,A)
print("\nN0 archived: median f_obs/f_pred by sample")
for sm in samples:
    m=samp==sm; print("  %-10s N=%2d  %.3f  (rms ln %.3f)"%(sm,m.sum(),np.exp(np.median(y[m]-np.log(f[m]))),np.std(y[m]-np.log(f[m]))))
print("\npredicted F(xi) archived: "+", ".join("xi=%.2f: f=%.3f"%(x,fF(np.array([x]),lam,n0,A)[0]) for x in (.03,.06,.1,.2,.3,.5)))
