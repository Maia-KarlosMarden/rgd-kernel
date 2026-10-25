"""RGD v6 candidate (reimagined radial sector).  NOT frozen for blind tests until the user decides; parameters below are the archived v5 values except p -> q=1.
Kernel:  K(x_eff) = 1 + lam exp(-x_eff)/(1-exp(-x_eff))^n,  n = n0 + A*xi.
Effective acceleration (geometric interpolation between the local and the halo acceleration):
    x_eff = x * (xi/y)^theta,   x = g_bar/g_dagger,  y = g_obs(r)/g_dagger,  xi = g_obs(r500)/g_dagger,  theta = q*Theta(xi),  Theta = [1+(xi_c/xi)^m]^-1.
theta=0 (galaxies): x_eff = x (local RAR kernel).  theta=1 (clusters, q=1): x_eff = f*xi -> f_bar independent of r inside a system.
Equivalent v5 form:  s = (g/g500)^theta  (v5 used s = (Delta/Delta500)^(p*Theta), p=0.5; locally p_eff = q(2-a)/(3-a), a = dlnM/dlnr).
Input: M [Msun] inside r [Mpc]; M500 [Msun], R500 [Mpc], z.  Output f_bar(<r)."""
import numpy as np
from scipy.optimize import brentq
G=6.674e-11; c=2.998e8; Msun=1.989e30; Mpc=3.0857e22
LAM,N0,A,Q,XIC,MW=0.999,0.5,0.4068,1.0,0.115,4
H=lambda z: 70e3/Mpc*np.sqrt(0.3*(1+z)**3+0.7)
def predict_fbar_v6(M,r_mpc,M500,R500_mpc,z,q=Q,xic=XIC,m=MW):
    gd=c*H(z)/(2*np.pi); r=r_mpc*Mpc; R5=R500_mpc*Mpc; g=G*M*Msun/r**2; g5=G*M500*Msun/R5**2
    xi=g5/gd; th=q/(1+(xic/xi)**m); n=N0+A*xi
    def h(f):
        x=f*g/gd; xe=x*(g5/g)**th
        return f*(1+LAM*np.exp(-xe)/(-np.expm1(-xe))**n)-1
    return brentq(h,1e-6,5.0),xi,th
if __name__=="__main__":
    import predict_frozen as pf
    M500,R500,z=8.3e14,1.38,0.078
    print("r/R500  M(<r)/M500  f_v6   f_v5(p=.5)")
    for w in (0.3,0.5,1.0,1.5):   # NFW c500~3 profile for the mass
        x=lambda w_: np.log(1+3*w_)-3*w_/(1+3*w_)
        Mr=M500*x(w)/x(1.0); f6=predict_fbar_v6(Mr,w*R500,M500,R500,z)[0]; f5=pf.predict_fbar(Mr,w*R500,M500,R500,z)[0]
        print("%.1f     %.3f       %.4f  %.4f"%(w,Mr/M500,f6,f5))
