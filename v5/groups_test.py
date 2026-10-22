import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.special import expit
from groups_data import LOV,SUN
G=6.674e-11; c=3e8; Msun=1.989e30; Mpc=3.0857e22; kpc=3.0857e19; lam=.999; A=.4068
def Hf(z,h,Om): return h*100e3/Mpc*np.sqrt(Om*(1+z)**3+1-Om)
grows=[]   # name,src,tier,z,M500,R500,M2500,R2500,fgas500,fgas2500 (Msun, m) ; H0 per source
for n,z,R5,M5,Mg5,R25,M25,Mg25 in LOV: grows.append(dict(name=n,src="L",tier=1,z=z,M5=M5*1e13*Msun,R5=R5*kpc,M25=M25*1e13*Msun,R25=R25*kpc,f5=Mg5*1e12/(M5*1e13),f25=Mg25*1e12/(M25*1e13),H=Hf(z,.7,.27)))
for n,z,t,R5,M5,f5,R25,f25 in SUN:
    H=Hf(z,.73,.24); rc=3*H**2/(8*np.pi*G); M25=2500*rc*4*np.pi/3*(R25*kpc)**3
    grows.append(dict(name=n,src="S",tier=t,z=z,M5=M5*1e13*Msun,R5=R5*kpc,M25=M25,R25=R25*kpc,f5=f5,f25=f25,H=H))
# unit sanity: M500 vs 500 rho_c (4pi/3) R500^3, and M2500 reported vs defined (Lovisari)
chk=[r['M5']/(500*3*r['H']**2/(8*np.pi*G)*4*np.pi/3*r['R5']**3) for r in grows]
chk25=[r['M25']/(2500*3*r['H']**2/(8*np.pi*G)*4*np.pi/3*r['R25']**3) for r in grows]
print("N=%d (Lovisari %d, Sun %d: tier1 %d, tier2 %d)"%(len(grows),len(LOV),len(SUN),sum(r['tier']==1 and r['src']=='S' for r in grows),sum(r['tier']==2 for r in grows)))
print("M500 / [500 rho_c 4pi/3 R500^3]: median %.3f range [%.2f,%.2f] ; M2500 same: median %.3f range [%.2f,%.2f]"%(np.median(chk),min(chk),max(chk),np.median(chk25),min(chk25),max(chk25)))
for r in grows:
    r['g500']=G*r['M5']/r['R5']**2; r['gd']=c*r['H']/(2*np.pi); r['xi']=r['g500']/r['gd']; r['w']=r['R25']/r['R5']
xi=np.array([r['xi'] for r in grows]); w=np.array([r['w'] for r in grows]); lM=np.log10([r['M5']/Msun for r in grows])
print("xi: %.3f-%.3f (median %.3f) ; logM500 %.2f-%.2f ; r2500/r500 median %.2f [%.2f,%.2f] ; fgas500 median %.3f ; fgas2500 median %.3f"%(xi.min(),xi.max(),np.median(xi),lM.min(),lM.max(),np.median(w),w.min(),w.max(),np.median([r['f5'] for r in grows]),np.median([r['f25'] for r in grows])))
def gresid(r,stars,sw,p=.5,at500=False):
    # sw: sat in [0,1]; evaluates log10(K*f_bar) at r2500 (s=5^(p*sw)) or at r500 (s=1)
    fb=r['f25']*(1+stars) if not at500 else r['f5']*(1+stars)
    M=r['M25'] if not at500 else r['M5']; R=r['R25'] if not at500 else r['R5']
    fbar=fb; x=fbar*G*M/R**2/r['gd']; s=1.0 if at500 else 5.0**(p*sw); n=.5+A*r['xi']
    xs=x/s; K=1+lam*np.exp(-xs)/(-np.expm1(-xs))**n
    return np.log10(K*fbar)
def stars_fn(mode,r):
    if mode=="paper": return 0.09
    if mode=="x2": return 0.09*2
    if mode=="scal": return min(0.39*(r['M5']/Msun/1e14)**-0.84,3.0)   # M*/Mgas scaling used for X-COP profile mode
    if mode=="none": return 0.0
rms=lambda a: np.sqrt(np.mean(np.square(a)))
print("\n=== A. kernel at R500 (s=1 by construction; checks n(xi) extrapolated to group xi) ===")
for mode in ("none","paper","x2","scal"):
    rr=np.array([gresid(r,stars_fn(mode,r),0,at500=True) for r in grows]); print(" stars=%-5s rms %.3f mean %+.3f | tier1+Lov only rms %.3f mean %+.3f | corr(resid,logM) %+.2f"%(mode,rms(rr),rr.mean(),rms(rr[[r['tier']==1 for r in grows]]),rr[[r['tier']==1 for r in grows]].mean(),np.corrcoef(rr,lM)[0,1]))
print("\n=== B. kernel at r2500 (Delta/Delta500=5): switch OFF (s=1) vs ON (s=5^0.5=2.24) ===")
for mode in ("none","paper","x2","scal"):
    off=np.array([gresid(r,stars_fn(mode,r),0) for r in grows]); on=np.array([gresid(r,stars_fn(mode,r),1) for r in grows])
    print(" stars=%-5s: OFF rms %.3f mean %+.3f | ON rms %.3f mean %+.3f"%(mode,rms(off),off.mean(),rms(on),on.mean()))
    for lo,hi in ((0,.15),(.15,.2),(.2,.25),(.25,.4)):
        mk=(xi>=lo)&(xi<hi)
        if mk.sum(): print("     xi [%.2f,%.2f) N=%2d: OFF mean %+.3f (rms %.3f) | ON mean %+.3f (rms %.3f)"%(lo,hi,mk.sum(),off[mk].mean(),rms(off[mk]),on[mk].mean(),rms(on[mk])))
np.save('groups_xi.npy',xi)
