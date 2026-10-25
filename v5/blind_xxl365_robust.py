# POST-HOC robustness (declared as such): concentration scale and a fixed-radius-independent check
import numpy as np, io, contextlib
import blind_xxl365_test as X
base=X.c200_DM14
for fac in (0.7,1.0,1.3):
    X.c200_DM14=lambda M,z,h=.7,fac=fac: fac*base(M,z,h)
    out=[]
    for d in X.rows:
        if not (np.isfinite(d['Mg']) and np.isfinite(d['r500']) and d['Mg']>0 and d['cls']==1): continue
        M5=X.M500_from_r(d['r500'],d['z']); r5=d['r500']*1e-3
        try: Mr=X.M_in(M5,r5,d['z'])
        except Exception: continue
        dd=dict(d,M5=M5,Mr=Mr,r5=r5); out.append(dd)
    for p in (0.5,0.8):
        lr=np.array([np.log10(X.pred(d,p=p)/(d['Mg']/d['Mr'])) for d in out if X.pred(d,p=p)>0])
        print("c x%.1f  p=%.1f  N=%d  median %+.3f dex (x%.2f)"%(fac,p,len(lr),np.median(lr),10**np.median(lr)))
