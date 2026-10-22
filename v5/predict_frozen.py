"""Frozen RGD v5 prescription (referee-proofed archive).  Apply to NEW samples before looking at the data.
Parameters are frozen here: lambda=0.999, A=0.4068, g_dagger(z)=c H(z)/(2 pi), H0=70, Om=0.3,
p=0.5, xi_c=0.115, m=4, stars: Kravtsov+2018 with phi=0.7 inside r2500 (1.0 inside r500 and r200).
Input (SI-free): M [Msun] enclosed mass (hydrostatic or lensing) inside r [Mpc]; M500 [Msun] and R500 [Mpc] of the system; z.
Output: predicted baryon fraction f_bar(<r) from the self-consistent equation f*K(f*g_obs/g_dagger)=1.
"""
import numpy as np
from scipy.optimize import brentq
G=6.674e-11; c=2.998e8; Msun=1.989e30; Mpc=3.0857e22
LAM,A,P,XIC,MW=0.999,0.4068,0.5,0.115,4
H=lambda z: 70e3/Mpc*np.sqrt(0.3*(1+z)**3+0.7)
def Mstar_tot(M500):  # Kravtsov+2018, Msun
    return 10**(12.71+0.59*(np.log10(M500)-14.5))
def predict_fbar(M,r_mpc,M500,R500_mpc,z,p=P,xic=XIC,m=MW):
    gd=c*H(z)/(2*np.pi); R5=R500_mpc*Mpc; r=r_mpc*Mpc; Mk=M*Msun; M5k=M500*Msun
    xi=G*M5k/R5**2/gd; theta=1/(1+(xic/xi)**m)
    dens=(Mk/M5k)*(R5/r)**3                      # Delta(r)/Delta_500
    s=dens**(p*theta); g=G*Mk/r**2; n=0.5+A*xi
    f=lambda fb: fb*(1+LAM*np.exp(-fb*g/gd/s)/(-np.expm1(-fb*g/gd/s))**n)-1
    return brentq(f,1e-6,5.0), xi, theta
if __name__=="__main__":
    # example: Abell 2029-like numbers, r = r2500-ish
    fb,xi,th=predict_fbar(M=3.0e14,r_mpc=0.6,M500=8.3e14,R500_mpc=1.38,z=0.078)
    print("f_bar(<r)=%.3f  xi=%.2f  Theta=%.2f"%(fb,xi,th))
