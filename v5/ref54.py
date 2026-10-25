"""PRE-REGISTERED (written before running).
Derivation: if s = g(r)/g500 (kernel argument radius-independent), then locally s = dd^p_eff with p_eff = dln g/dln dd = (2-a)/(3-a), a = dlnM(<r)/dlnr.
Test (zero free parameters for the formula): per-system p_i (ref51 machinery, 254+ systems) vs p_pred,i = (2-a_i)/(3-a_i) with a_i from an NFW profile (Dutton&Maccio14 concentration,
solved so that M(<R500)=M500) evaluated at the system's measurement radius.  Models: constant p; p=p_pred (0 params); p=b0+b1*p_pred (2 params); each with hierarchical tau.
Criteria: AIC and leave-one-sample-out score as in ref51.  Caveat declared: for XXL, M(<500 kpc) itself came from the same NFW (so a_i is self-consistent but p_i depends on the gas data only).
Also print p_pred at r2500, r500, r200 for c=3.5."""
import io,contextlib,sys,numpy as np,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import ref51 as Q
R=Q.R; X=Q.X; D=R.D; Mpc=R.Mpc; Ms=R.Ms
from scipy.optimize import brentq
nfw=lambda x: np.log(1+x)-x/(1+x)
def a_of(x): return x**2/((1+x)**2*nfw(x))        # dlnM/dlnr for NFW, x=c r/R200
print("NFW c200=3.5 (R500=0.65*R200 approx): effective p=(2-a)/(3-a)")
for lab,xx in (("r2500",3.5*0.30),("r500",3.5*0.65),("r200",3.5*1.0)):
    a=a_of(xx); print("   %-6s x=%.2f a=%.3f p_eff=%.3f"%(lab,xx,a,(2-a)/(3-a)))
rows=Q.rows
name2d={d['name']:d for d in D}
def p_pred(d):
    M5=d['M5']/Ms; R5=d['R5']/Mpc; z=d['z']; r=d['r']/Mpc
    rho=3*X.Hz2(z)**2/(8*np.pi*X.G)
    f=lambda lM: 10**lM*X.nfw(X.c200_DM14(10**lM,z)*R5/((3*10**lM*Ms/(4*np.pi*200*rho))**(1/3)/Mpc))/X.nfw(X.c200_DM14(10**lM,z))-M5
    lM=brentq(f,np.log10(M5)+0.01,np.log10(M5)+1.0); M200=10**lM; R200=(3*M200*Ms/(4*np.pi*200*rho))**(1/3)/Mpc; c=X.c200_DM14(M200,z)
    a=a_of(c*r/R200); return (2-a)/(3-a),a,c
pp=[];keep=[]
for i,rw in enumerate(rows):
    d=name2d.get(rw[0])
    if d is None: continue
    try: v=p_pred(d)
    except Exception: continue
    pp.append(v[0]); keep.append(i)
keep=np.array(keep); pp=np.array(pp)
y=Q.y[keep]; s=Q.s[keep]; samp=Q.samp[keep]
print("\nN=%d  p_pred range %.2f..%.2f  median %.2f ; observed p_i median %.2f"%(len(y),pp.min(),pp.max(),np.median(pp),np.median(y)))
for sm in sorted(set(samp)):
    m=samp==sm; print("  %-18s N=%3d  p_i median %.2f   p_pred median %.2f"%(sm,m.sum(),np.median(y[m]),np.median(pp[m])))
def nll(res,ss,tau): v=ss**2+tau**2; return 0.5*np.sum(np.log(2*np.pi*v)+res**2/v)
from scipy.optimize import minimize_scalar
def tau_fit(res,ss): o=minimize_scalar(lambda lt: nll(res,ss,np.exp(lt)),bounds=(-6,1),method='bounded'); return np.exp(o.x),o.fun
out={}
# constant
def fitm(mode,idx):
    if mode=='const':
        w=1/(s[idx]**2+0.02); mu=np.sum(w*y[idx])/np.sum(w)
        for _ in range(30):
            t,_=tau_fit(y[idx]-mu,s[idx]); w=1/(s[idx]**2+t**2); mu=np.sum(w*y[idx])/np.sum(w)
        return (mu,),t
    if mode=='formula': t,_=tau_fit(y[idx]-pp[idx],s[idx]); return (),t
    if mode=='lin':
        X_=np.column_stack([np.ones(idx.sum()),pp[idx]]); t=0.15
        for _ in range(30):
            w=1/(s[idx]**2+t**2); A=X_.T@(X_*w[:,None]); b=np.linalg.solve(A,X_.T@(w*y[idx])); t,_=tau_fit(y[idx]-X_@b,s[idx])
        return tuple(b),t
def mu_of(mode,th,idx):
    if mode=='const': return np.full(idx.sum(),th[0])
    if mode=='formula': return pp[idx]
    return th[0]+th[1]*pp[idx]
allm=np.ones(len(y),bool)
print("\nmodel      k  tau    nll     AIC    LOSO     params")
for mode,k in (('const',1),('formula',0),('lin',2)):
    th,t=fitm(mode,allm); v=nll(y-mu_of(mode,th,allm),s,t); aic=2*v+2*(k+1)
    lo=0
    for sm in sorted(set(samp)):
        te=samp==sm; tr=~te; tht,tt=fitm(mode,tr); lo-=nll(y[te]-mu_of(mode,tht,te),s[te],tt)
    print("%-9s %d  %.3f  %6.1f  %6.1f  %7.1f   %s"%(mode,k,t,v,aic,lo,np.round(th,3)))
