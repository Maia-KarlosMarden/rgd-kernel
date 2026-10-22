import numpy as np
from scipy.optimize import brentq, least_squares
from mantz_data import MANTZ
G=6.674e-11; c=3e8; Msun=1.989e30; kpc=3.0857e19; Mpc=3.0857e22; lam=.999; A=.4068
mf=lambda x: np.log(1+x)-x/(1+x)
Hz=lambda z: 70e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
def M500_from_M2500(M25,r25,z,c500):
    rc=3*Hz(z)**2/(8*np.pi*G)
    def f(lR5):
        R5=10**lR5; return M25*mf(c500*r25/R5)/mf(c500)*0+ (lambda M5: M5*mf(c500*r25/R5)/mf(c500)-M25)(500*rc*4*np.pi/3*R5**3)
    R5=10**brentq(f,np.log10(r25*1.2),np.log10(r25*6)); return 500*rc*4*np.pi/3*R5**3,R5
def fpred(M25,r25,z,xi,gd,p,s_on=True):
    g=G*M25/r25**2; s=5.0**p if s_on else 1.0; n=.5+A*xi
    K=lambda f: 1+lam*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**n
    return brentq(lambda f: f*K(f)-1,1e-4,5.0)
D=[(n,z,R*kpc,M*1e14*Msun,f,e) for n,z,R,M,f,e in MANTZ]
z=np.array([d[1] for d in D]); fo=np.array([d[4] for d in D]); ef=np.array([d[5] for d in D]); M25=np.array([d[3] for d in D])
def run(c500,p,stars,verbose=False):
    xi=[];pr=[];prF=[]
    for n,zz,r25,m25,f,e in D:
        M5,R5=M500_from_M2500(m25,r25,zz,c500); x=G*M5/R5**2/(c*Hz(zz)/(2*np.pi)); xi.append(x)
        pr.append(fpred(m25,r25,zz,x,c*Hz(zz)/(2*np.pi),p)); prF.append(fpred(m25,r25,zz,x,c*Hz(0)/(2*np.pi),p))
    return np.array(xi),np.array(pr),np.array(prF)
def chi2_free_norm(model,obs,err):
    # obs = N*model ; weighted least squares in ln
    lo=np.log(obs); le=err/obs; lm=np.log(model); w=1/le**2; lnN=np.sum(w*(lo-lm))/np.sum(w); r=(lo-lm-lnN)/le; return np.sum(r**2),np.exp(lnN)
def slope(y,x,err):
    w=1/err**2; X=np.vstack([np.ones_like(x),x]).T; cov=np.linalg.inv(X.T@(w[:,None]*X)); b=cov@X.T@(w*y); return b[1],np.sqrt(cov[1,1])
print("40 clusters, z=0.078-1.063. Observed: f_gas in the 0.8-1.2 r2500 shell (not cumulative) -> normalisation left free; we test the z-SHAPE.")
for c500,p,stars in ((3.0,.5,.09),(3.0,.8,.09),(2.5,.5,.09),(4.0,.5,.09),(3.0,.5,.0)):
    xi,pr,prF=run(c500,p,stars)
    fb_obs=fo*(1+stars)
    ln_r=np.log(fb_obs/pr)
    c2R,NR=chi2_free_norm(pr,fb_obs,ef*(1+stars)); c2F,NF=chi2_free_norm(prF,fb_obs,ef*(1+stars)); c2C,NC=chi2_free_norm(np.ones_like(pr),fb_obs,ef*(1+stars))
    bz,ez=slope(np.log(fb_obs),np.log(1+z),ef/fo)
    print("\n c500=%.1f p=%.1f stars=%.2f | xi %.2f-%.2f | RGD predicted f_bar(<r2500) %.3f-%.3f (median %.3f; z<0.3 median %.3f, z>0.5 median %.3f)"%(c500,p,stars,xi.min(),xi.max(),pr.min(),pr.max(),np.median(pr),np.median(pr[z<.3]),np.median(pr[z>.5])))
    print("  observed f_bar median %.3f (z<0.3 %.3f ; z>0.5 %.3f) ; observed dln f/dln(1+z) = %+.2f +- %.2f"%(np.median(fb_obs),np.median(fb_obs[z<.3]),np.median(fb_obs[z>.5]),bz,ez))
    bz_p,_=slope(np.log(pr),np.log(1+z),ef/fo); bz_F,_=slope(np.log(prF),np.log(1+z),ef/fo)
    print("  predicted dln f/dln(1+z): RGD g_dag(z) = %+.2f ; g_dag frozen = %+.2f ; LCDM-like = 0"%(bz_p,bz_F))
    print("  chi2 (free norm, 39 dof): RGD g_dag(z) %.1f (N=%.2f) | g_dag frozen %.1f (N=%.2f) | constant (LCDM-like) %.1f (N=%.2f)"%(c2R,NR,c2F,NF,c2C,NC))
    print("  absolute check, no normalisation: median obs/pred = %.2f  [shell f_gas exceeds cumulative f_gas, so >1 expected]"%np.exp(np.median(ln_r)))

print("\n=== mass dependence (c500=3, p=0.5, stars 0.09) ===")
xi,pr,prF=run(3.0,.5,.09); fb=fo*1.09; lM=np.log10(M25/Msun)
w=1/(ef*1.09/fb)**2
def wslope(y,x): X=np.vstack([np.ones_like(x),x]).T; cov=np.linalg.inv(X.T@(w[:,None]*X)); b=cov@X.T@(w*y); return b[1],np.sqrt(cov[1,1])
print(" dln f_bar / dlog10 M2500: observed %+.2f +- %.2f ; RGD predicted %+.2f"%(*wslope(np.log(fb),lM),wslope(np.log(pr),lM)[0]))
print(" corr(ln obs, ln pred) = %.2f ; Spearman-like check by M2500 tercile (obs median / pred median):"%np.corrcoef(np.log(fb),np.log(pr))[0,1])
q=np.quantile(lM,[1/3,2/3])
for lo,hi,nm in ((-9,q[0],"low M"),(q[0],q[1],"mid M"),(q[1],9,"high M")):
    mk=(lM>=lo)&(lM<hi); print("   %-6s N=%2d  obs %.3f  pred %.3f"%(nm,mk.sum(),np.median(fb[mk]),np.median(pr[mk])))
print(" z-trend at fixed mass: use high-mass half only (log M2500 > median):")
mk=lM>np.median(lM); zz=z[mk]
b1=np.polyfit(np.log(1+zz),np.log(fb[mk]),1)[0]; b2=np.polyfit(np.log(1+zz),np.log(pr[mk]),1)[0]; b3=np.polyfit(np.log(1+zz),np.log(prF[mk]),1)[0]
print("   dln f/dln(1+z): obs %+.2f ; RGD %+.2f ; frozen g_dag %+.2f"%(b1,b2,b3))
