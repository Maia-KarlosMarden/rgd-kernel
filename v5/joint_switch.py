import numpy as np, warnings; warnings.filterwarnings("ignore")
from scipy.special import expit
from scipy.optimize import least_squares
cd={}; exec(open('kernel_scenarios.py').read().split("# ---------- check NFW R500")[0],cd)
src=open('sparc_bridge.py').read().split("base=resid(")[0]
exec(src)
Dc=cd['build'](radii_mpc=[.5,1,1.5]); Gc=cd['G']; Ac=cd['A']; lamc=cd['lam']
xc=Dc['fb']*Dc['g']/Dc['gd']; densc=(Dc['g']/(Gc*Dc['M5']/Dc['R5']**2))/Dc['w']
xi_c=Dc['xi']; M5c=Dc['M5']/1.989e30 if Dc['M5'].max()>1e30 else Dc['M5']
print("clusters: N=%d  xi range %.3f-%.3f ; log M500 %.2f-%.2f ; Delta/Delta500 median %.1f [%.1f,%.1f]"%(xi_c.size,xi_c.min(),xi_c.max(),np.log10(M5c.min()),np.log10(M5c.max()),np.median(densc),densc.min(),densc.max()))
def Kc(s): xs=xc/s; return 1+lamc*np.exp(-xs)/(-np.expm1(-xs))**(0.5+Ac*xi_c)
def rc(sat,p): return np.log10(Kc(np.clip(densc,1e-6,None)**(p*sat))*Dc['fb'])
print(" cluster rms: no radial term (sat=0) %.4f ; full (sat=1,p=.5) %.4f ; p fitted sat=1: "%(np.sqrt(np.mean(rc(0,.5)**2)),np.sqrt(np.mean(rc(1,.5)**2))),end="")
fp=least_squares(lambda q: rc(1,q[0]).ravel(),[.5]); print("p=%.3f rms %.4f"%(fp.x[0],np.sqrt(np.mean(fp.fun**2))))
mg=inset
def rs(sat,p,idx=None): 
    s=np.clip(dens,1e-6,None)**(p*sat); return (np.log10(gb*K(gb/gd,s=s))-np.log10(gobs))[mg]
def m2ll(rcl,rsp): return rcl.size*np.log(np.mean(rcl**2))+rsp.size*np.log(np.mean(rsp**2))
base=m2ll(rc(0,.5),rs(np.zeros(len(R)),.5)); allon=m2ll(rc(1,.5),rs(np.ones(len(R)),.5))
print("-2lnL: all-off (s=1 everywhere) %.0f ; all-on (p=.5 everywhere) %.0f"%(base,allon))
Mhalo_pt={}
for model in ("A","B"):
    xi_h=prep(model); 
    Mh_pt=np.zeros(len(R))
    for g in gl: Mh_pt[ids==g]=halo(Mb[g],model)
    print("\n=== SPARC halo model %s: xi_halo %.4f-%.3f ; logMh %.1f-%.1f ==="%(model,xi_h[mg].min(),xi_h[mg].max(),np.log10(Mh_pt[mg].min()),np.log10(Mh_pt[mg].max())))
    for var,vc,vs in (("xi_halo",xi_c,xi_h),("M500",M5c,Mh_pt)):
        print(" variable =",var)
        lvc=np.log(vc); lvs=np.log(np.where(vs>0,vs,1))
        grid=np.linspace(min(lvc.min(),lvs[mg].min())-0.5,max(lvc.max(),lvs[mg].max())+0.5,400)
        for mm in (1,2,4,8,32,200):
            best=(1e99,)
            for l0 in grid:
                v=m2ll(rc(expit(mm*(lvc-l0)),.5),rs(expit(mm*(lvs-l0)),.5))
                if v<best[0]: best=(v,l0)
            l0=best[1]; satc=expit(mm*(lvc-l0)); sats=expit(mm*(lvs-l0))
            rcl=rc(satc,.5); rsp=rs(sats,.5)
            print("  m=%3d: -2lnL=%7.0f (dAIC vs all-off %+.0f) | thr %s=%.3g | rms clusters %.4f SPARC %.4f | frac sat>0.5: clusters %.2f SPARC %.2f"%(mm,best[0],best[0]+4-base,var,np.exp(l0),np.sqrt(np.mean(rcl**2)),np.sqrt(np.mean(rsp**2)),np.mean(satc>.5),np.mean(sats[mg]>.5)))
    # best sharp switch on xi_h: free p
    print(" free p at m=32 on xi_halo:")
    lvc=np.log(xi_c); lvs=np.log(np.where(xi_h>0,xi_h,1))
    f=lambda q: np.concatenate([rc(expit(32*(lvc-q[1])),q[0]).ravel()/1.0, rs(expit(32*(lvs-q[1])),q[0])])
    gridl=np.linspace(np.log(0.05),np.log(0.3),60); bb=(1e99,)
    for l0 in gridl:
        for p in (0.3,0.4,0.5,0.6,0.8):
            v=m2ll(rc(expit(32*(lvc-l0)),p),rs(expit(32*(lvs-l0)),p))
            if v<bb[0]: bb=(v,l0,p)
    print("   best -2lnL %.0f thr xi=%.3f p=%.2f"%(bb[0],np.exp(bb[1]),bb[2]))
