import numpy as np
from scipy.optimize import brentq
G=6.674e-11; c=3e8; Msun=1.989e30; Mpc=3.0857e22; lam=.999; A=.4068
mf=lambda x: np.log(1+x)-x/(1+x)
H=lambda z,const=False: 70e3/Mpc*(1 if const else np.sqrt(.3*(1+z)**3+.7))
def r_delta(M500,z,D,c500=3.0):
    rc=3*H(z)**2/(8*np.pi*G); R5=(3*M500*Msun/(4*np.pi*500*rc))**(1/3); 
    f=lambda lr: M500*Msun*mf(c500*10**lr/R5)/mf(c500)-D*rc*4*np.pi/3*(10**lr)**3
    r=10**brentq(f,np.log10(R5*.05),np.log10(R5*3)); return r,M500*Msun*mf(c500*r/R5)/mf(c500),R5
def fbar_pred(M500,z,D,gdz,xi_for_n,s=1.0):
    r,M,R5=r_delta(M500,z,D); g=G*M/r**2
    K=lambda f: 1+lam*np.exp(-f*g/gdz/s)/(-np.expm1(-f*g/gdz/s))**(.5+A*xi_for_n)
    sol=brentq(lambda f: f*K(f)-1,1e-3,1.0); return sol
print("RGD-predicted baryon fraction within r2500 (Delta=2500; switch ON, s=5^0.5) for a fixed M500 and different z; LCDM expects f_gas(z)~const.")
print(" Compare g_dagger(z)=cH(z)/2pi (RGD) with g_dagger frozen at z=0 (control)")
for M in (3e14,1e15):
    print(" M500=%.0e:"%M)
    for z in (0,.3,.6,1.0,1.5):
        rc=3*H(z)**2/(8*np.pi*G); R5=(3*M*Msun/(4*np.pi*500*rc))**(1/3); xi=G*M*Msun/R5**2/(c*H(z)/(2*np.pi))
        a=fbar_pred(M,z,2500,c*H(z)/(2*np.pi),xi,s=5**.5); b=fbar_pred(M,z,2500,c*H(0)/(2*np.pi),xi,s=5**.5)
        print("   z=%.1f  xi=%.2f  f_bar(<r2500): g_dag(z)=%.4f ; g_dag frozen=%.4f ; ratio to z=0 (g_dag(z)) = %.3f"%(z,xi,a,b,a/fbar_pred(M,0,2500,c*H(0)/(2*np.pi),G*M*Msun/((3*M*Msun/(4*np.pi*500*3*H(0)**2/(8*np.pi*G)))**(2/3))/(c*H(0)/(2*np.pi)),s=5**.5)))
print("\nBTFR/RAR zero-point: V^4 ∝ g_dagger*M_b -> at fixed M_b, dlogV = 0.25 log10 E(z)")
for z in (.5,1,1.5,2,2.5): print("  z=%.1f E=%.2f  dlogV=%+.3f dex (%.0f%%)"%(z,np.sqrt(.3*(1+z)**3+.7),.25*np.log10(np.sqrt(.3*(1+z)**3+.7)),100*(np.sqrt(.3*(1+z)**3+.7)**.25-1)))
print("\nThreshold mass M500 where xi crosses xi_c (z-dependence of the switch)")
for z in (0,.5,1,1.5,2):
    out=[]
    for xc in (.09,.14):
        Ms=np.logspace(11,15,4000); v=[]
        for Mx in Ms:
            rc=3*H(z)**2/(8*np.pi*G); R5=(3*Mx*Msun/(4*np.pi*500*rc))**(1/3); v.append(G*Mx*Msun/R5**2/(c*H(z)/(2*np.pi)))
        out.append("%.1e"%Ms[np.argmin(abs(np.array(v)-xc))])
    print("  z=%.1f  xi_c=0.09 -> %s ; xi_c=0.14 -> %s Msun"%(z,out[0],out[1]))
