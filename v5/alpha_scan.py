import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.optimize import brentq
from allen_data import ALLEN
from mantz_data import MANTZ
G=6.674e-11; c=3e8; Ms=1.989e30; kpc=3.0857e19; Mpc=3.0857e22; lam=.999; A=.4068
mf=lambda x: np.log(1+x)-x/(1+x)
Hz=lambda z: 70e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
gd=lambda z,al: c*Hz(0)*(Hz(z)/Hz(0))**al/(2*np.pi)
def R5_of(M25,r25,z,c500):
    rc=3*Hz(z)**2/(8*np.pi*G)
    f=lambda l: 500*rc*4*np.pi/3*(10**l)**3*mf(c500*r25/10**l)/mf(c500)-M25
    R5=10**brentq(f,np.log10(r25*1.2),np.log10(r25*6)); return 500*rc*4*np.pi/3*R5**3,R5
def fpred(m,r,xi,g_d,p):
    g=G*m/r**2; s=5.0**p
    h=lambda f: f*(1+lam*np.exp(-f*g/g_d/s)/(-np.expm1(-f*g/g_d/s))**(.5+A*xi))-1
    try: return brentq(h,1e-9,5.0)
    except ValueError: return np.nan
def chi2(model,obs,err):
    if not np.all(np.isfinite(model)): return np.inf,np.nan
    lo=np.log(obs); le=err/obs; w=1/le**2; lnN=np.sum(w*(lo-np.log(model)))/np.sum(w); return np.sum(((lo-np.log(model)-lnN)/le)**2),np.exp(lnN)
# datasets: (z, r2500, M2500, f_obs, err)
rc_=lambda z: 3*Hz(z)**2/(8*np.pi*G)
AL=[(a[1],a[2]*kpc,2500*rc_(a[1])*4*np.pi/3*(a[2]*kpc)**3,a[3],a[4]) for a in ALLEN]
MA=[(m[1],m[2]*kpc,m[3]*1e14*Ms,m[4],m[5]) for m in MANTZ]
ST=.09; C500=3.0
def scan(data,als,p=.5):
    out=[]
    for al in als:
        pr=[];fo=[];eo=[]
        for z,r,m,f,e in data:
            M5,R5=R5_of(m,r,z,C500); xi=G*M5/R5**2/gd(z,al); pr.append(fpred(m,r,xi,gd(z,al),p)); fo.append(f*(1+ST)); eo.append(e*(1+ST))
        out.append(chi2(np.array(pr),np.array(fo),np.array(eo)))
    return out
als=np.round(np.arange(0,1.51,.25),2)
for p in (.5,.46):
    ra=scan(AL,als,p); rm=scan(MA,als,p)
    print("p=%.2f   alpha:  "%p+"  ".join("%6.2f"%a for a in als))
    print("   Allen chi2 (41 dof):  "+"  ".join("%6.1f"%x[0] for x in ra)); print("   Mantz chi2 (39 dof):  "+"  ".join("%6.1f"%x[0] for x in rm))
    tot=np.array([a[0]+b[0] for a,b in zip(ra,rm)]); i=tot.argmin(); print("   Sum (Allen+Mantz; overlap!):"+"  ".join("%6.1f"%x for x in tot)+"  -> min at alpha=%.2f"%als[i])
# fine scan for 1-sigma on Allen alone and Mantz alone, p=0.5
fine=np.round(np.arange(0,1.51,.05),2)
for nm,d in (("Allen",AL),("Mantz",MA)):
    r=np.array([x[0] for x in scan(d,fine)]); j=r.argmin(); ok=fine[r<r[j]+1]; ok4=fine[r<r[j]+4]
    print("%s (p=0.5): best alpha=%.2f chi2=%.1f ; 1-sigma [%.2f, %.2f] ; 2-sigma [%.2f, %.2f] ; chi2(alpha=1)=%.1f chi2(alpha=0)=%.1f"%(nm,fine[j],r[j],ok.min(),ok.max(),ok4.min(),ok4.max(),r[list(fine).index(1.0)],r[list(fine).index(0.0)]))
