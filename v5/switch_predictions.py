import numpy as np
G=6.674e-11; c=3e8; Msun=1.989e30; Mpc=3.0857e22; lam=.999; A=.4068
def H(z): return 70e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
def R500(M,z): rc=3*H(z)**2/(8*np.pi*G); return (3*M*Msun/(4*np.pi*500*rc))**(1/3)
def xi(M,z): return G*M*Msun/R500(M,z)**2/(c*H(z)/(2*np.pi))
mf=lambda x: np.log(1+x)-x/(1+x)
def effect(M,z,w=0.5,c500=2.7):
    x_=xi(M,z); R5=R500(M,z); g500=G*M*Msun/R5**2; gd=c*H(z)/(2*np.pi)
    mr=mf(c500*w)/mf(c500); g=g500*mr/w**2; fb=1.09*0.0708*(M/1e13)**0.135
    dens=(g/g500)/w; xx=fb*g/gd; n=.5+A*x_
    K=lambda s:1+lam*np.exp(-xx/s)/(-np.expm1(-xx/s))**n
    return np.log10(K(dens**0.5)/K(1.0))        # change in predicted mass if term ON vs OFF
print("xi = g500/g_dagger(z) for M500 and z (xi ~ (M*H)^(1/3)); regime under switch xi_c in [0.09,0.34]:  ON if xi>0.34, OFF if xi<0.09, '?' between")
print("%8s"%"M500", *["z=%.1f        "%z for z in (0,0.5,1,2)])
for M in (3e12,1e13,3e13,1e14,3e14,1e15):
    row=[]
    for z in (0,.5,1,2):
        x_=xi(M,z); tag="ON " if x_>.34 else ("OFF" if x_<.09 else "?  ")
        row.append("%.2f %s %+.2f"%(x_,tag,effect(M,z)))
    print("%8.0e"%M,"  ".join(row))
print("(last number = dex difference in predicted M(<0.5 R500) between switch ON and OFF; measurable if >~0.05 dex)")
print("\nMass at which xi crosses xi_c, by redshift:")
for z in (0,.5,1,2):
    out=[]
    for xc in (.09,.16,.34):
        Ms=np.logspace(11,16,2000); v=np.array([xi(M,z) for M in Ms]); out.append("%.1e"%Ms[np.argmin(abs(v-xc))])
    print(" z=%.1f  xi_c=0.09/0.16/0.34 -> M500c = %s Msun"%(z,", ".join(out)))
print("\nXLSSC 122 (M500=1.66e14, z=1.98): xi=%.3f -> ON ; effect at w=0.5: %+.2f dex"%(xi(1.66e14,1.98),effect(1.66e14,1.98)))
