import numpy as np, sys
from blind_xxl_efeds import *
# 1) eFEDS vs HSC-XXL observed gas, same mass/z
ef=[(3.0,.138,.022),(5.6,.211,.043),(8.4,.263,.052),(13.5,.310,.065),(13.3,.470,.049),(22.5,.697,.053),(25.1,.483,.063),(27.7,.318,.070)]
print("obs eFEDS f_gas / HSC-XXL f_gas (same M500, z=0.3):")
for m,z,fg in ef:
    g,s=xxl_frac(m*1e13,.3); print("  M=%5.1fe13  eFEDS %.3f  XXL %.3f  ratio %.2f   RGD_gas/eFEDS %.2f"%(m,fg,g,fg/g,(fbar(m*1e13,z)-xxl_frac(m*1e13,z)[1])/fg))
# 2) stellar fraction: Kravtsov vs HSC-XXL
print("f_* Kravtsov(frozen) vs HSC-XXL measured:")
for M in [1.5e13,3e13,1e14,3e14,1e15]:
    print("  M500=%.1e  Kravtsov %.4f  HSC-XXL %.4f"%(M,pf.Mstar_tot(M)/M,xxl_frac(M,.3)[1]))
# 3) best p at r200 vs Popesso eq5 for M200 in 3e13..3e14, and at r500 pred vs eq4 (non-indep.)
h=.6774; z=.15
def c200(M):
    a=.520+(.905-.520)*np.exp(-.617*z**1.21); b=-.101+.026*z
    return 10**(a+b*np.log10(M/(1e12/h)))
nfw=lambda x: np.log(1+x)-x/(1+x)
def setup(M200):
    R200=rdelta(M200,z,200,67.74,.3089); c=c200(M200)
    r5=brentq(lambda r: nfw(c*r/R200)/nfw(c)*M200-2.5*(r/R200)**3*M200,0.2*R200,R200)
    return R200,r5,M200*nfw(c*r5/R200)/nfw(c)
Ms=np.array([3e13,5e13,1e14,2e14,3e14])
def resid(p,xic=.115,m=4):
    out=[]
    for M200 in Ms:
        R200,R500,M500=setup(M200)
        pg=pf.predict_fbar(M200,R200,M500,R500,z,p=p,xic=xic,m=m)[0]-pf.Mstar_tot(M500)/M200
        out.append(np.log10(pg/(2.09e-6*M200**.33)))
    return np.array(out)
ps=np.arange(.3,1.01,.05); rms=[np.sqrt(np.mean(resid(p)**2)) for p in ps]
i=int(np.argmin(rms)); print("Popesso r200 (3e13-3e14): best p=%.2f rms=%.3f dex ; p=0.5 rms=%.3f"%(ps[i],rms[i],np.sqrt(np.mean(resid(.5)**2))))
print(" resid(dex) p=0.5:",np.round(resid(.5),3)," p=%.2f:"%ps[i],np.round(resid(ps[i]),3))
# r500 Popesso eq4 over same masses
def r500res(M200):
    R200,R500,M500=setup(M200); pg=pf.predict_fbar(M500,R500,M500,R500,z)[0]-pf.Mstar_tot(M500)/M500
    return np.log10(pg/(2.23e-7*M500**.39))
print("Popesso r500 resid (dex), M200=3e13..3e14:",np.round([r500res(M) for M in Ms],3))
# with measured stars (HSC-XXL) instead of Kravtsov at r500
def r500res2(M200):
    R200,R500,M500=setup(M200); s=xxl_frac(M500,z)[1]; pg=pf.predict_fbar(M500,R500,M500,R500,z)[0]-s
    return np.log10(pg/(2.23e-7*M500**.39))
print("  same with HSC-XXL stars:",np.round([r500res2(M) for M in Ms],3))
