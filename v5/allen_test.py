import numpy as np
from scipy.optimize import brentq
from allen_data import ALLEN
G=6.674e-11; c=3e8; Msun=1.989e30; kpc=3.0857e19; Mpc=3.0857e22; lam=.999; A=.4068
mf=lambda x: np.log(1+x)-x/(1+x)
Hz=lambda z: 70e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
print("N=%d clusters, z=%.3f-%.3f (cumulative f_gas within r2500)"%(len(ALLEN),ALLEN[0][1],ALLEN[-1][1]))
z=np.array([a[1] for a in ALLEN]); r25=np.array([a[2] for a in ALLEN])*kpc; fo=np.array([a[3] for a in ALLEN]); ef=np.array([a[4] for a in ALLEN])
M25=np.array([2500*3*Hz(zz)**2/(8*np.pi*G)*4*np.pi/3*r**3 for zz,r in zip(z,r25)])
def M500(c500):
    xi=[]
    for zz,r,m in zip(z,r25,M25):
        rc=3*Hz(zz)**2/(8*np.pi*G); f=lambda l: 500*rc*4*np.pi/3*(10**l)**3*mf(c500*r/10**l)/mf(c500)-m
        R5=10**brentq(f,np.log10(r*1.2),np.log10(r*6)); M5=500*rc*4*np.pi/3*R5**3; xi.append(G*M5/R5**2/(c*Hz(zz)/(2*np.pi)))
    return np.array(xi)
def fsolve(m,r,xi,gd,p):
    g=G*m/r**2; s=5.0**p; K=lambda f: 1+lam*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**(.5+A*xi); return brentq(lambda f: f*K(f)-1,1e-4,5.0)
def wsl(y,x,e):
    w=1/e**2; X=np.vstack([np.ones_like(x),x]).T; cv=np.linalg.inv(X.T@(w[:,None]*X)); b=cv@X.T@(w*y); return b[1],np.sqrt(cv[1,1])
def chi2fn(model,obs,err):
    lo=np.log(obs); le=err/obs; w=1/le**2; lnN=np.sum(w*(lo-np.log(model)))/np.sum(w); return np.sum(((lo-np.log(model)-lnN)/le)**2),np.exp(lnN)
lM=np.log10(M25/Msun); print("M2500 = %.2e-%.2e Msun ; median %.2e"%(M25.min()/Msun,M25.max()/Msun,np.median(M25)/Msun))
for c500 in (3.0,2.5,4.0):
    xi=M500_=M500(c500)
    for p in (.5,.46,.8):
        for st in (.09,.0):
            pr=np.array([fsolve(m,r,x,c*Hz(zz)/(2*np.pi),p) for m,r,x,zz in zip(M25,r25,xi,z)]); prF=np.array([fsolve(m,r,x,c*Hz(0)/(2*np.pi),p) for m,r,x,zz in zip(M25,r25,xi,z)])
            fb=fo*(1+st); e=ef*(1+st)
            if not (c500==3.0 or (p==.5 and st==.09)): continue
            if c500!=3.0 and p!=.5: continue
            c_R,N_R=chi2fn(pr,fb,e); c_F,N_F=chi2fn(prF,fb,e); c_C,N_C=chi2fn(np.ones_like(pr),fb,e)
            bM,eM=wsl(np.log(fb),lM,e/fb); bMp=np.polyfit(lM,np.log(pr),1)[0]; bz,ez=wsl(np.log(fb),np.log(1+z),e/fb); bzp=np.polyfit(np.log(1+z),np.log(pr),1)[0]; bzF=np.polyfit(np.log(1+z),np.log(prF),1)[0]
            print("\nc500=%.1f p=%.2f stars=%.2f | xi %.2f-%.2f | pred f_bar %.3f-%.3f (median %.3f) obs median %.3f -> obs/pred (median) %.2f"%(c500,p,st,xi.min(),xi.max(),pr.min(),pr.max(),np.median(pr),np.median(fb),np.median(fb/pr)))
            print("  chi2 (free norm, 41 dof): RGD g_dag(z) %.1f (N=%.2f) | g_dag frozen %.1f | constant %.1f (N=%.3f)"%(c_R,N_R,c_F,c_C,N_C))
            print("  mass slope dln f/dlog M2500: obs %+.2f+-%.2f ; RGD %+.2f ; z slope dln f/dln(1+z): obs %+.2f+-%.2f ; RGD %+.2f ; frozen %+.2f"%(bM,eM,bMp,bz,ez,bzp,bzF))

print("\n=== which ingredient generates the mass slope? (c500=3, stars=0.09, p=0.5) ===")
xi=M500(3.0); fb=fo*1.09; e=ef*1.09
def pred(Aeff,p,nfix=None):
    out=[]
    for m,r,x,zz in zip(M25,r25,xi,z):
        gd=c*Hz(zz)/(2*np.pi); g=G*m/r**2; s=5.0**p; n=(.5+Aeff*x) if nfix is None else nfix
        K=lambda f: 1+lam*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**n; out.append(brentq(lambda f: f*K(f)-1,1e-4,5.0))
    return np.array(out)
for Aeff in (0.4068,0.3,0.2,0.1,0.0,-0.1):
    pr=pred(Aeff,.5); c2,N=chi2fn(pr,fb,e); print(" n=0.5+%.2f xi: pred slope %+.2f (obs -0.04+-0.09) ; chi2 %.1f (free norm N=%.2f) ; median obs/pred %.2f"%(Aeff,np.polyfit(lM,np.log(pr),1)[0],c2,N,np.median(fb/pr)))
print(" fixed n(xi)=0.5 (galactic value):",end=" ")
for p in (.5,.6,.7):
    pr=pred(0,p,nfix=.5); c2,N=chi2fn(pr,fb,e); print("p=%.1f slope %+.2f chi2 %.1f obs/pred %.2f;"%(p,np.polyfit(lM,np.log(pr),1)[0],c2,np.median(fb/pr)),end=" ")
print()
