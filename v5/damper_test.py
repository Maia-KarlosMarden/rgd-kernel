import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
Om=.3;Ob=.049;Or=9.2e-5;OL=1-Om-Or;lam=.999;A=.4068;Kc=Om/Ob
E2=lambda a: Om*a**-3+Or*a**-4+OL
dlnE=lambda a:(-1.5*Om*a**-3-2*Or*a**-4)/E2(a)
Omb=lambda a: Ob*a**-3/E2(a)
c=2.998e5;H0=70.;RH=lambda a:c/(H0*np.sqrt(E2(a)))
expit=lambda u: 1/(1+np.exp(-np.clip(u,-50,50)))
def K(x,n=.5,s=1.): x=max(x,1e-30); return 1+lam*np.exp(-x/s)/(-np.expm1(-x/s))**n
def Keff(x,Db,Dv=30.,m=4.):          # Db = baryon overdensity rho_b/rho_c ; D=1 bound, 0 cosmic
    D=expit(m*np.log(max(Db,1e-30)/Dv)); return Kc*(1-D)+K(x)*D
def run(Rc,d0=1e-3,mode="damper"):
    a0=1/1101.
    def f(t,y):
        a=np.exp(t); d,dp=y
        if mode=="lcdm": src=1.5*Om*a**-3/E2(a)*d
        else:
            x=np.pi*Omb(a)*abs(d)*a*Rc/RH(a); Db=Omb(a)*(1+d); src=1.5*Omb(a)*Keff(x,Db)*d
        return [dp,-(2+dlnE(a))*dp+src]
    s=solve_ivp(f,[np.log(a0),0],[d0,d0],rtol=1e-9,atol=1e-14,method="LSODA",dense_output=True); return s
print("1) Linear growth with damper D(rho_b/rho_c; Dv=30, m=4) and cosmic plateau Kc=Om/Ob=%.2f"%Kc)
for Rc in (1,10,100):
    s=run(Rc); l=run(Rc,mode="lcdm")
    print("   R_c=%3g Mpc: delta(z=0)  damper %.4g  LCDM %.4g   max |dlnD/dlna - LCDM| = %.1e"%(Rc,s.y[0,-1],l.y[0,-1],max(abs(s.sol(t)[1]/s.sol(t)[0]-l.sol(t)[1]/l.sol(t)[0]) for t in np.linspace(np.log(1/1101.),0,200))))
print("   -> linear growth equals LCDM by construction while rho_b/rho_c < Dv (here at z=0 delta_b=1 gives rho_b/rho_c = %.3f << 30):"%(Ob*2))
print("      the switch only engages after collapse (rho_b >~ 30 rho_c, i.e. total Delta ~ 200), as intended.")
print("\n2) Non-trivial consequence: cumulative f_bar at r200 in the Vikhlinin clusters (RGD, p=0.5) with and without the cap 1/Kc=%.3f"%(1/Kc))
G=6.674e-11; ccc=3e8; Ms=1.989e30; kpc=3.0857e19; Mpc=3.0857e22
H=lambda z:71e3/Mpc*np.sqrt(.3*(1+z)**3+.7); mf=lambda x:np.log(1+x)-x/(1+x)
# name, z, r500 kpc, c500, M500 (1e14)
V=[("A133",.0569,998,3.31,3.14),("A383",.1883,956,4.36,3.10),("A478",.0881,1359,3.75,7.83),("A907",.1603,1117,3.56,4.71),("A1413",.1429,1339,3.06,7.78),("A1795",.0622,1283,3.73,6.57),("A1991",.0592,753,4.66,1.28),("A2029",.0779,1380,4.13,8.29),("A2390",.2302,1448,1.65,10.88),("MKW4",.0199,635,2.69,.74)]
print("   %-6s %5s | f_bar(r500) | f_bar(r200) uncapped  capped | r200/r500"%("name","xi"))
tot=[];
for n,z,r5,c5,M5 in V:
    rc=3*H(z)**2/(8*np.pi*G); R5=r5*kpc; M5k=M5*1e14*Ms
    # NFW through (R5,M5): find r200 with M(<r)=200 rc 4pi/3 r^3
    rs=R5/c5; Mn=lambda r: M5k*mf(r/rs)/mf(c5)
    R2=brentq(lambda r: Mn(r)-200*rc*4*np.pi/3*r**3,R5,3*R5); M2=Mn(R2)
    gd=ccc*H(z)/(2*np.pi); xi=G*M5k/R5**2/gd
    def f_of(m,r,Dd):
        s=Dd**.5; g=G*m/r**2
        return brentq(lambda f: f*(1+lam*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**(.5+A*xi))-1,1e-4,5)
    f5=f_of(M5k,R5,1.0); f2=f_of(M2,R2,0.4)    # Delta/Delta500 = 1 at r500, 0.4 at r200
    f2c=max(f2,1/Kc); tot.append((f5,f2,f2c))
    print("   %-6s %5.2f |    %.3f    |    %.3f              %.3f   |  %.2f"%(n,xi,f5,f2,f2c,R2/R5))
t=np.array(tot); print("   median: r500 %.3f ; r200 uncapped %.3f ; r200 capped %.3f ; cosmic Ob/Om = %.3f"%(np.median(t[:,0]),np.median(t[:,1]),np.median(t[:,2]),1/Kc))
