import numpy as np, collections
from scipy.optimize import least_squares
import os
f=os.environ.get('SPARC_TABLE2','SPARC_Table2.txt')
rows=[l.split() for l in open(f).read().splitlines()[25:] if l.strip()]
ids=np.array([r[0] for r in rows]); a=np.array([[float(x) for x in r[1:]] for r in rows])
D,R,V,eV,Vg,Vd,Vb=a[:,:7].T
print("points",len(R),"galaxies",len(set(ids)))
kpc=3.0857e19; gd=3e8*70e3/3.0857e22/(2*np.pi)
def gbar(Yd,Yb): return (Vg*np.abs(Vg)+Yd*Vd**2+Yb*Vb**2)*1e6/(R*kpc)
gobs=V**2*1e6/(R*kpc)
ok=(V>0)&(R>0)&(eV/V<0.1)
def K(x,s=1.0,n=0.5,lam=0.999):
    xs=x/s; return 1+lam*np.exp(-xs)/(-np.expm1(-xs))**n
def stats(pred,m):
    r=np.log10(pred[m])-np.log10(gobs[m]); return np.sqrt(np.mean(r**2)),r.mean()
gb=gbar(0.5,0.7); m=ok&(gb>0)
print("N used",m.sum(),"| g_dagger(z=0)=%.3e"%gd)
print("RGD galactic n=1/2 (s=1): rms %.3f dex, mean %+.3f"%stats(gb*K(gb/gd),m))
mond=gb/(1-np.exp(-np.sqrt(gb/1.2e-10))); print("McGaugh RAR (a0=1.2e-10): rms %.3f mean %+.3f"%stats(mond,m))
print("Newtonian: rms %.3f mean %+.3f"%stats(gb,m))
# exploratory radial term: mean-density ratio relative to outermost point of each galaxy
ref=np.zeros_like(R); 
for g in set(ids):
    i=np.where(ids==g)[0]; j=i[np.argmax(R[i])]; ref[i]=(V[j]**2/R[j]**2)
dens=np.where(ref>0,(V**2/R**2)/np.where(ref>0,ref,1),1)
for p in (0.25,0.5,1.0):
    print("exploratory s=(Delta/Delta_out)^p, p=%.2f: rms %.3f mean %+.3f"%((p,)+stats(gb*K(gb/gd,s=dens**p),m)))
fit=least_squares(lambda q: np.log10((gb*K(gb/gd,s=dens**q[0]))[m])-np.log10(gobs[m]),[0.1]); print("fitted p=%.3f rms %.3f"%(fit.x[0],np.sqrt(np.mean(fit.fun**2))))
fy=least_squares(lambda q: np.log10((gbar(q[0],0.7)*0+gbar(q[0],0.7)*K(gbar(q[0],0.7)/gd))[m])-np.log10(gobs[m]),[0.5]); print("s=1 with global fitted Y_disk=%.2f rms %.3f"%(fy.x[0],np.sqrt(np.mean(fy.fun**2))))
# by regime
lg=np.log10(gb)
for lo,hi in ((-13,-11.5),(-11.5,-10.5),(-10.5,-9)):
    mm=m&(lg>=lo)&(lg<hi); print(" log gbar [%g,%g): n=%d rgd resid mean %+.3f rms %.3f"%(lo,hi,mm.sum(),*stats(gb*K(gb/gd),mm)[::-1]))
# dwarf vs massive
vm=collections.defaultdict(float)
for g in set(ids): i=ids==g; vm[g]=V[i&ok].max() if (i&ok).any() else 0
big=np.array([vm[g]>150 for g in ids]); small=np.array([vm[g]<70 for g in ids])
for nm,mk in (("Vmax>150",big),("Vmax<70",small)): print(nm,"n=%d rms %.3f mean %+.3f"%((m&mk).sum(),*stats(gb*K(gb/gd),m&mk)))
