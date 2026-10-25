"""FROZEN 2026-10-25 -- competing predictions for blind tests.  DO NOT EDIT after the freeze tag (v6-freeze-20261025); new versions go in new files.
Calibration source: ref58.py (systematics-aware fit, sd_delta=0.20, 331 points, 6 samples).  Everything else is the archived v5 prescription.
v6:  x_eff = x * (xi/y)^(q*Theta(xi)),  q = 1.26   (equivalently s = (g/g500)^(q Theta))
v5:  s = (Delta/Delta500)^(p*Theta(xi)),  p = 0.61
Common: K = 1 + lam exp(-x/s)/(1-exp(-x/s))^n, lam=.999, n=.5+.4068 xi, g_dag=cH(z)/2pi (H0=70, Om=.3), Theta=[1+(xi_c/xi)^m]^-1, xi_c=.115, m=4,
 stars: Kravtsov+2018 total stellar mass, phi=0.7 inside r2500 and 1.0 inside r500 and r200 (applied to the observed baryon fraction by the user of this module).
Amplitude systematics attached to every prediction: ln-offset per sample ~ N(0, 0.20^2); intrinsic scatter tau: v6 0.151, v5 0.198 (ln f).
Input: M [Msun] inside r [Mpc]; M500 [Msun], R500 [Mpc], z.  Output f_bar(<r)."""
import numpy as np
from scipy.optimize import brentq
G=6.674e-11; c=2.998e8; Msun=1.989e30; Mpc=3.0857e22
FROZEN=dict(lam=0.999,n0=0.5,A=0.4068,xic=0.115,m=4,q_v6=1.26,p_v5=0.61,sd_offset=0.20,tau_v6=0.151,tau_v5=0.198,phi_r2500=0.7,phi_r500_r200=1.0)
H=lambda z: 70e3/Mpc*np.sqrt(0.3*(1+z)**3+0.7)
def Mstar_tot(M500):   # Kravtsov+2018, Msun
    return 10**(12.71+0.59*(np.log10(M500)-14.5))
def _solve(M,r_mpc,M500,R500_mpc,z,law):
    F=FROZEN; gd=c*H(z)/(2*np.pi); r=r_mpc*Mpc; R5=R500_mpc*Mpc; Mk=M*Msun; M5k=M500*Msun
    g=G*Mk/r**2; g5=G*M5k/R5**2; xi=g5/gd; th=1/(1+(F['xic']/xi)**F['m']); n=F['n0']+F['A']*xi
    s=(g/g5)**(F['q_v6']*th) if law=='v6' else ((Mk/M5k)*(R5/r)**3)**(F['p_v5']*th)
    h=lambda f: f*(1+F['lam']*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**n)-1
    return brentq(h,1e-6,5.0)
def predict_v6(M,r_mpc,M500,R500_mpc,z): return _solve(M,r_mpc,M500,R500_mpc,z,'v6')
def predict_v5(M,r_mpc,M500,R500_mpc,z): return _solve(M,r_mpc,M500,R500_mpc,z,'v5')
if __name__=="__main__":
    print("v6/v5 f_bar for a 8.3e14 Msun cluster, z=0.078, NFW c500~3:")
    x=lambda w: np.log(1+3*w)-3*w/(1+3*w)
    for w in (0.3,0.5,1.0,1.5):
        Mr=8.3e14*x(w)/x(1.0); print("  r/R500=%.1f  v6 %.4f   v5 %.4f"%(w,predict_v6(Mr,w*1.38,8.3e14,1.38,0.078),predict_v5(Mr,w*1.38,8.3e14,1.38,0.078)))
