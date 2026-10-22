import io, contextlib
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('groups_test.py').read())
xi_g=np.array([r['xi'] for r in grows]); N_g=len(grows)
with contextlib.redirect_stdout(buf): exec(open('joint_switch.py').read().split("Mhalo_pt={}")[0])
from scipy.optimize import least_squares
def rg(satv,p,mode):
    return np.array([gresid(r,stars_fn(mode,r),satv[i],p=p) for i,r in enumerate(grows)])
print("GROUP SAMPLE: N=%d, xi %.3f-%.3f ; every one has Delta(r2500)/Delta500 = 5 exactly"%(N_g,xi_g.min(),xi_g.max()))
print("\n(1) free p at r2500 (sat=1) by xi bin  [p=0 -> switch OFF; clusters need p~0.45-0.5]")
for mode in ("paper","x2","scal"):
    s_=" stars=%-5s:"%mode
    for lo,hi in ((.13,.18),(.18,.22),(.22,.30),(.13,.30)):
        mk=np.where((xi_g>=lo)&(xi_g<hi))[0]
        f=least_squares(lambda q: np.array([gresid(grows[i],stars_fn(mode,grows[i]),1,p=q[0]) for i in mk]),[.3],bounds=(-.5,2))
        s_+="  xi[%.2f,%.2f) N=%2d p=%.2f (rms %.3f);"%(lo,hi,len(mk),f.x[0],rms(f.fun))
    print(s_)
print("\n(2) joint 3-dataset scan of the sharp switch (m=32), p=1/2: -2lnL(xi_c) = Nc ln rms_c^2 + Ns ln rms_s^2 + Ng ln rms_g^2")
grid=np.exp(np.linspace(np.log(0.04),np.log(0.5),220)); lg=np.log(xi_g)
for model in ("A","B"):
    xi_h=prep(model); lvs=np.log(np.where(xi_h>0,xi_h,1)); lvc=np.log(xi_c)
    for mode in ("paper","x2","scal"):
        tot=[];parts=[]
        for xcut in grid:
            lx=np.log(xcut); rcl=rc(expit(32*(lvc-lx)),.5); rsp=rs(expit(32*(lvs-lx)),.5); rgp=rg(expit(32*(lg-lx)),.5,mode)
            tot.append(rcl.size*np.log(np.mean(rcl**2))+rsp.size*np.log(np.mean(rsp**2))+rgp.size*np.log(np.mean(rgp**2))); parts.append((rms(rcl),rms(rsp),rms(rgp)))
        tot=np.array(tot); i=np.argmin(tot); ok_=grid[tot<tot[i]+4]
        offg=rg(np.zeros(N_g),.5,mode); ong=rg(np.ones(N_g),.5,mode)
        print(" SPARC halo %s, stars=%-5s: best xi_c=%.3f ; 2-sigma-ish interval [%.3f, %.3f] ; rms at best: clusters %.3f SPARC %.3f groups %.3f | groups all-OFF %.3f all-ON %.3f"%(model,mode,grid[i],ok_.min(),ok_.max(),*parts[i],rms(offg),rms(ong)))
print("\n(3) What each dataset alone allows (stars=paper, SPARC halo A):")
xi_h=prep("A"); lvs=np.log(np.where(xi_h>0,xi_h,1)); lvc=np.log(xi_c)
for nm,fn in (("SPARC",lambda lx: rs(expit(32*(lvs-lx)),.5)),("groups",lambda lx: rg(expit(32*(lg-lx)),.5,"paper"))):
    v=np.array([fn(np.log(x)).size*np.log(np.mean(fn(np.log(x))**2)) for x in grid]); ok_=grid[v<v.min()+4]; print(" %-6s alone: xi_c in [%.3f, %.3f]  (rms at xi_c=0.10: %.3f; 0.20: %.3f; 0.30: %.3f)"%(nm,ok_.min(),ok_.max(),rms(fn(np.log(.1))),rms(fn(np.log(.2))),rms(fn(np.log(.3)))))
