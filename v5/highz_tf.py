import numpy as np
from scipy.special import i0,i1,k0,k1
G=4.30091e-6          # kpc (km/s)^2 / Msun
c=2.998e5; Mpc_kpc=1e3; lam=.999
H0=70.                # km/s/Mpc
E=lambda z: np.sqrt(.3*(1+z)**3+.7)
conv=1e6/3.0857e19    # (km/s)^2/kpc -> m/s^2
gdag=lambda z: c*1e3*(H0*1e3/3.0857e22)*E(z)/(2*np.pi)       # m/s^2
def K(x,n=.5): 
    x=np.maximum(x,1e-30); return 1+lam*np.exp(-x)/(-np.expm1(-x))**n
def vcirc_max(Mb,Re,z,n=.5,frozen=False):
    Rd=Re/1.678; r=np.linspace(.2,12,600)*Rd; y=r/(2*Rd)
    v2b=2*G*Mb/Rd*y**2*(i0(y)*k0(y)-i1(y)*k1(y))          # Freeman: 4 pi G S0 Rd y^2(...), S0=Mb/(2 pi Rd^2)
    gb=v2b/r*conv; x=gb/(gdag(0 if frozen else z)); v2=K(x,n)*v2b
    j=np.argmax(v2); return np.sqrt(v2[j]), x[j]
OBS=[("z~0.9",0.9,10.62,4.8,239),("z~1.5",1.5,10.97,4.9,263),("z~2.3",2.3,10.89,4.0,260)]
print("RGD circular velocity of the KMOS3D median galaxies (exp. disk, all baryons, n=1/2, s=1)")
print("%-6s %5s | logMb  Re  | V_obs | V_RGD(g†(z))  V_RGD(g† frozen) | V_RGD(z=0 same M,R) | pred/obs [z-dep, frozen] | x at peak"%("bin","E(z)"))
rows=[]
for nm,z,lMb,Re,Vo in OBS:
    Mb=10**lMb; v1,x1=vcirc_max(Mb,Re,z); v2,_=vcirc_max(Mb,Re,z,frozen=True); v0,_=vcirc_max(Mb,Re,0.0)
    rows.append((nm,z,Vo,v1,v2,v0)); print("%-6s %5.2f | %.2f %.1f | %5d |   %6.1f         %6.1f           |    %6.1f          |   %.2f   %.2f   | %.2f"%(nm,E(z),lMb,Re,Vo,v1,v2,v0,v1/Vo,v2/Vo,x1))
print("\nlog10 shift of V_pred relative to z=0 (same M_b,R_e): " + "; ".join("%s %+.3f dex"%(r[0],np.log10(r[3]/r[5])) for r in rows))
print("bTFR-style zero-point shift  dlogM = 3.75*dlogV (Lelli slope): "+"; ".join("%s pred %+.2f"%(r[0],-3.75*np.log10(r[3]/r[5])) for r in rows))
print("Observed zero-point vs local (Ubler Table/§3.3): z~0.9 -0.44 ; z~2.3 -0.27 (dex in logM at fixed V)")
for n in (.5,.55): print("sensitivity n=%.2f:"%n, ["%s %.1f"%(nm,vcirc_max(10**lMb,Re,z,n=n)[0]) for nm,z,lMb,Re,Vo in OBS])
print("dependence on an assumed gas distribution (Re x1.3): ", ["%s %.1f"%(nm,vcirc_max(10**lMb,1.3*Re,z)[0]) for nm,z,lMb,Re,Vo in OBS])
