# Referee point 3.4: is a sharp xi-step needed? Free logistic width m, threshold xi_c and p; smooth power-law alternative; equal-parameter comparison
import io, contextlib, numpy as np, warnings; warnings.filterwarnings("ignore")
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('groups_test.py').read())
with contextlib.redirect_stdout(buf):
    exec(open('stars_test.py').read().split('print("\\nALLEN')[0])
xi_g=np.array([r['xi'] for r in grows]); N_g=len(grows); lg=np.log(xi_g)
with contextlib.redirect_stdout(buf): exec(open('joint_switch.py').read().split("Mhalo_pt={}")[0])
from scipy.special import expit
KRA=np.array([.7*Mstar_tot(r['M5']/Msun)*Msun/(r['f25']*r['M25']) for r in grows])
def rg(sat,p,mode):
    st=KRA if mode=="kra" else np.array([stars_fn(mode,r) for r in grows])
    return np.array([gresid(r,st[i],sat[i],p=p) for i,r in enumerate(grows)])
R_=lambda a: float(np.sqrt(np.mean(np.square(a))))
import sys
NE=(13,141,43) if len(sys.argv)<2 else (39,2803,43)   # effective independent units: clusters, galaxies, groups (default); points if argv given
def total(satc,sats,satg,p,mode):
    rcl=rc(satc,p); rsp=rs(sats,p); rgp=rg(satg,p,mode)
    return NE[0]*np.log(np.mean(rcl**2))+NE[1]*np.log(np.mean(rsp**2))+NE[2]*np.log(np.mean(rgp**2)),(R_(rcl),R_(rsp),R_(rgp))
lvc=np.log(xi_c)
XC=np.exp(np.linspace(np.log(0.04),np.log(0.6),70)); PP=np.arange(0.3,1.01,0.05); MM=[1,2,3,4,6,8,16,32,200]
print("Datasets: X-COP 39 pts (xi %.2f-%.2f), groups 43 (xi %.2f-%.2f), SPARC (xi_halo imputed from the halo model, <0.1)"%(xi_c.min(),xi_c.max(),xi_g.min(),xi_g.max()))
print('weights N_eff (clusters,SPARC,groups) =',NE)
for model in ("A","B"):
    xi_h=prep(model); lvs=np.log(np.where(xi_h>0,xi_h,1))
    for mode in ("kra","paper"):
        print("\n=== SPARC halo %s, group stars=%s ==="%(model,mode))
        allOff=total(np.zeros(xi_c.shape),np.zeros(len(R)),np.zeros(N_g),0.5,mode); allOn=total(np.ones(xi_c.shape),np.ones(len(R)),np.ones(N_g),0.5,mode)
        print(" reference (p=.5): switch OFF everywhere %.0f (rms c/s/g %.3f/%.3f/%.3f) ; ON everywhere %.0f (%.3f/%.3f/%.3f)"%(allOff[0],*allOff[1],allOn[0],*allOn[1]))
        best_m={}
        for m in MM:
            b=(1e99,)
            for xcut in XC:
                lx=np.log(xcut); sc=expit(m*(lvc-lx)); ss=expit(m*(lvs-lx)); sg=expit(m*(lg-lx))
                for p in PP:
                    v,parts=total(sc,ss,sg,p,mode)
                    if v<b[0]: b=(v,xcut,p,parts)
            best_m[m]=b
        mb=min(best_m,key=lambda k:best_m[k][0])
        for m in MM:
            v,xcb,p,parts=best_m[m]; print("  m=%3d: -2lnL %8.1f (d vs best %+6.1f) xi_c=%.3f p=%.2f | rms clusters %.3f SPARC %.3f groups %.3f"%(m,v,v-best_m[mb][0],xcb,p,*parts))
        # smooth power-law p(xi)=p0*(xi/0.4)^q , no threshold, 2 parameters
        b=(1e99,)
        for q in np.arange(0,1.51,0.05):
            for p0 in np.arange(0.3,1.31,0.05):
                v,parts=total((xi_c/0.4)**q,(np.where(xi_h>0,xi_h,1e-9)/0.4)**q,(xi_g/0.4)**q,p0,mode)
                if v<b[0]: b=(v,p0,q,parts)
        print("  POWER-LAW p(xi)=p0 (xi/0.4)^q (2 par, no threshold): -2lnL %.1f (d vs best step/logistic %+.1f) p0=%.2f q=%.2f | rms c/s/g %.3f/%.3f/%.3f"%(b[0],b[0]-best_m[mb][0],*b[1:3],*b[3]))
        v32=best_m[32][0]; print("  => parameters: sharp step (xi_c,p)=2 ; logistic (xi_c,m,p)=3 ; power-law (p0,q)=2.  AIC: step %.1f  logistic(best m=%d) %.1f  power-law %.1f"%(v32+4,mb,best_m[mb][0]+6,b[0]+4))
