import numpy as np
from scipy.optimize import least_squares, brentq
from scipy.stats import spearmanr
exec(open('sparc_mbar.py').read().split("base=np.sqrt")[0])
Mpc_=3.0857e22; H0=70e3/Mpc_; rhoc=3*H0**2/(8*np.pi*G)   # kg/m3
def halo(Mb_,model):
    if model=="A": return Mb_/0.10                                   # const f_bar(R500)=0.10
    f=lambda lMh: np.log10(Mb_)-(np.log10(0.04)+1.4*lMh-0.4*12)-0*lMh   # Mb = 0.04*Mh*(Mh/1e12)^0.4
    return 10**brentq(lambda l: f(l)*1.0, 8, 16)
gal_idx={g:i for i,g in enumerate(gl)}
def prep(model):
    xi_h=np.zeros(len(R)); 
    for g in gl:
        Mh=halo(Mb[g],model)*Msun; R5=(3*Mh/(4*np.pi*500*rhoc))**(1/3); xi_h[ids==g]=(G*Mh/R5**2)/gd
    return xi_h
Mtot=V**2*1e6*(R*kpc)/G/Msun            # total enclosed mass from rotation curve (Msun)
dens=(3*Mtot*Msun/(4*np.pi*(R*kpc)**3))/(500*rhoc)   # Delta(r)/Delta500 (rho_bar(<r)/500 rho_c)
inset=np.isin(ids,gl)&m
def resid(xi_h,p,X0,mask): 
    sat=1-np.exp(-xi_h/X0); s=np.clip(dens,1e-6,None)**(p*sat); return np.log10(gb*K(gb/gd,s=s))[mask]-np.log10(gobs[mask])
rms_=lambda r: np.sqrt(np.mean(r**2))
base=resid(np.zeros(len(R)),0,0.052,inset); print("baseline s=1: rms %.3f (N=%d)"%(rms_(base),inset.sum()))
print("Delta(r)/Delta500 over SPARC points: median %.0f, 5-95%% [%.0f, %.0f] (X-COP clusters: ~1-10)"%(np.median(dens[inset]),*np.percentile(dens[inset],[5,95])))
for model in ("A","B"):
    xi_h=prep(model); xg=np.array([xi_h[ids==g][0] for g in gl]); lMh=np.log10([halo(Mb[g],model) for g in gl])
    print("\n=== halo model %s : log Mh %.1f-%.1f ; xi_halo median %.4f range [%.4f, %.3f]; sat(0.052) median %.2f ==="%(model,lMh.min(),lMh.max(),np.median(xg),xg.min(),xg.max(),np.median(1-np.exp(-xg/0.052))))
    print(" (1) paper prediction p=1/2, X0=0.052, nothing fitted: rms %.3f"%rms_(resid(xi_h,0.5,0.052,inset)))
    f=least_squares(lambda q: resid(xi_h,0.5,10**q[0],inset),[-3.0],bounds=(-8,1)); X0=10**f.x[0]
    print(" (2) p=1/2, X0 fitted: log10 X0=%.2f (X0=%.2e vs 0.052) rms %.3f ; sat at median xi=%.2e"%(f.x[0],X0,rms_(f.fun),np.median(1-np.exp(-xg/X0))))
    # profile in X0
    print("     rms vs X0:",", ".join("%.0e:%.3f"%(x,rms_(resid(xi_h,0.5,x,inset))) for x in (1e-1,1e-2,1e-3,1e-4,1e-5,1e-6)))
    # CV
    rng=np.random.default_rng(2); fo=np.array([rng.integers(0,5) if False else 0 for _ in ids]); fold={g:rng.integers(0,5) for g in gl}; fo=np.array([fold.get(g,-1) for g in ids]); cv=[]
    for k in range(5):
        tr=inset&(fo!=k); te=inset&(fo==k); q=least_squares(lambda q: resid(xi_h,0.5,10**q[0],tr),[-3.0],bounds=(-8,1)).x; cv+=list(resid(xi_h,0.5,10**q[0],te))
    print(" (3) 5-fold CV rms with fitted X0: %.3f (baseline %.3f)"%(rms_(np.array(cv)),rms_(base)))
    # free p, sat=1: does any mass bin want p>0?
    print(" (4) free p (sat=1) by halo-mass bin:")
    for lo,hi in ((7,10.5),(10.5,11.5),(11.5,12.2),(12.2,15)):
        mk=np.zeros(len(R),bool)
        for g,l in zip(gl,lMh):
            if lo<=l<hi: mk|=(ids==g)
        mk&=inset
        if mk.sum()<30: continue
        r0=rms_(resid(xi_h,0,0.052,mk)); ff=least_squares(lambda q: resid(np.full(len(R),1e9),q[0],0.052,mk),[0.05],bounds=(-1,2))
        print("     logMh [%.1f,%.1f): N=%d  p_best=%.3f  rms(p=0)=%.3f rms(p_best)=%.3f  rms(p=0.5)=%.3f"%(lo,hi,mk.sum(),ff.x[0],r0,rms_(ff.fun),rms_(resid(np.full(len(R),1e9),0.5,0.052,mk))))
