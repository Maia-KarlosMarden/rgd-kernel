"""PRE-REGISTERED (written before running).
Question: can luminosity selection inflate the mean f_gas of a sample enough to explain RGD(>)eFEDS?
Model: at fixed true M500, (ln Lx, ln Mg) bivariate normal; sigma_L=0.73, sigma_g=0.39, cov(L,g)=0.28 (HSC-XXL Table 2, Table 1),
  r=0.98 (central) and r=0.80 (low end of the quoted covariance error). Selection: Lx > cut, detection fraction P_det at fixed M.
  Bias on mean ln Mg = r*sigma_g*lambda(alpha), lambda = phi(alpha)/(1-Phi(alpha)), P_det = 1-Phi(alpha).
  Stars are ignored in the bias (corr(Lx,M*)~0.18, small); bias applied to the gas fraction only (log-mean; mean of f_gas itself gets +sigma_g^2/2 on both sides, cancels in the ratio).
Outputs: (a) mean f_gas boost vs P_det; (b) the P_det that would reconcile each eFEDS bin with RGD (assuming eFEDS is the unbiased population mean and RGD is the mean of a selected sample).
No survey-specific cut is assumed.
"""
import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq
sg=0.39
def boost(P,r):
    a=norm.isf(P); lam=norm.pdf(a)/P
    return np.exp(r*sg*lam)
print("P_det   boost(r=0.98)  boost(r=0.80)")
for P in [1,.9,.7,.5,.3,.2,.1,.05,.02,.01]:
    print("%5.2f   %.2f          %.2f"%(P,boost(P,.98),boost(P,.80)))
ef=[(3.0,.138,2.58),(5.6,.211,1.62),(8.4,.263,1.51),(13.5,.310,1.36),(13.3,.470,1.83),(22.5,.697,1.96),(25.1,.483,1.64),(27.7,.318,1.48)]
print("\nRGD_gas/eFEDS_gas ratio per bin and P_det that reconciles it")
for m,z,R in ef:
    out=[]
    for r in (.98,.80):
        try: out.append(brentq(lambda P: boost(P,r)-R,1e-6,.9999))
        except ValueError: out.append(float('nan'))
    print("  M=%5.1fe13  ratio %.2f   P_det(r=.98)=%.3f   P_det(r=.80)=%.4f"%(m,R,out[0],out[1]))
# versus XXL: RGD/obs_f_b ~0.97-1.06 -> boost needed ~1 -> P_det ~1
print("\nRGD vs HSC-XXL (ratio 0.97-1.06): reconciling boost ~1 => P_det ~ 1 (no selection excess needed).")
