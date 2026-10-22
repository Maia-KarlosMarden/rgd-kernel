exec(open('ref34.py').read().split("lvc=np.log(xi_c)")[0])
import sys
for mode in ("kra","paper"):
    for p in (0.45,0.5,0.8):
        r=rg(np.ones(N_g),p,mode); print("groups stars=%s p=%.2f: rms %.3f mean %+.3f"%(mode,p,R_(r),r.mean()))
r=rg(np.ones(N_g),0.0,"kra"); print("groups kra p=0 (OFF): rms %.3f mean %+.3f"%(R_(r),r.mean()))
from scipy.optimize import minimize_scalar
for mode in ("kra","paper"):
    o=minimize_scalar(lambda p: R_(rg(np.ones(N_g),p,mode)),bounds=(0,1.5),method="bounded"); print("best p",mode,o.x,o.fun)
# bootstrap CI of best p for groups (kra)
rng=np.random.default_rng(3)
for mode in ("kra","paper"):
    bp=[]
    for _ in range(300):
        ix=rng.integers(0,N_g,N_g)
        f=lambda p: R_(np.array([gresid(grows[i],(KRA[i] if mode=="kra" else stars_fn("paper",grows[i])),1,p=p) for i in ix]))
        bp.append(minimize_scalar(f,bounds=(0,1.5),method="bounded").x)
    print(mode,"group p 95%% CI [%.2f, %.2f]"%tuple(np.percentile(bp,[2.5,97.5])))
