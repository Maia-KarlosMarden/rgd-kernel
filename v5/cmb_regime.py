import numpy as np
from scipy.integrate import solve_ivp
c=2.998e5; H0=70.; Om=.3; Ob=.049; Or=9.2e-5; OL=1-Om-Or; lam=.999
E=lambda a: np.sqrt(Om*a**-3+Or*a**-4+OL)
dlnE=lambda a: (-1.5*Om*a**-3-2*Or*a**-4)/E(a)**2
Omb=lambda a: Ob*a**-3/E(a)**2
Omm=lambda a: Om*a**-3/E(a)**2
RH=lambda a: c/(H0*E(a))                    # Mpc physical
def K(x,n=.5,s=1.): 
    x=max(x,1e-30); return 1+lam*np.exp(-x/s)/(-np.expm1(-x/s))**n
# x = g_bar/g_dagger for a spherical overdensity of baryons: g=(1/2)Omb H^2 d R ; g_dag=cH/2pi  -> x=pi*Omb*d*R/RH
def x_of(a,d,Rc): return np.pi*Omb(a)*abs(d)*(a*Rc)/RH(a)
print("1) Regime map: x = g_bar/g_dagger for a linear baryon perturbation (comoving R)")
print("   and K(x) (n=1/2, s=1) versus the boost a CDM universe needs (Om/Ob = %.1f)"%(Om/Ob))
for z,d in [(1100,1e-4),(1100,1e-3),(100,1e-3),(10,1e-2),(1,.1),(0,.1),(0,1.)]:
    a=1/(1+z); print("  z=%5d delta_b=%.0e: "%(z,d)+"  ".join("R_c=%4g Mpc: x=%.1e K=%.0f"%(R,x_of(a,d,R),K(x_of(a,d,R))) for R in (0.1,1,10,100)))
print("\n2) Toy growth, baryons only (Ob), K(x(delta)) vs LCDM (Om, K=1), same initial delta at z=1100")
def run(d0,Rc,mode):
    a0=1/1101.; 
    def f(t,y):
        a=np.exp(t); d,dp=y
        src=1.5*(Omm(a) if mode=="lcdm" else Omb(a)*K(x_of(a,d,Rc)))*d
        return [dp, -(2+dlnE(a))*dp+src]
    # growing-mode initial condition d' = d (radiation era approx; irrelevant for ratios)
    s=solve_ivp(f,[np.log(a0),0],[d0,d0],rtol=1e-8,atol=1e-14,method="LSODA"); return s.y[0,-1]
for Rc in (1,10,100):
    for d0 in (1e-5,1e-4,1e-3):
        dl=run(d0,Rc,"lcdm"); dr=run(d0,Rc,"rgd")
        print("  R_c=%4g Mpc  d_i=%.0e : LCDM d(z=0)=%.3g   RGD(baryon-only) d(z=0)=%.3g   ratio=%.2g"%(Rc,d0,dl,dr,dr/dl))
print("\n3) Linearity test: LCDM d_f/d_i is amplitude independent; for RGD:")
for Rc in (10,):
    r=[run(d0,Rc,"rgd")/d0 for d0 in (1e-6,1e-5,1e-4,1e-3)]; print("   R_c=%g Mpc: d_f/d_i for d_i=1e-6..1e-3 :"%Rc,["%.3g"%v for v in r], "  LCDM:", "%.3g"%(run(1e-5,Rc,'lcdm')/1e-5))
