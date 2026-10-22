# Referee point 3.2: sensitivity of the radial-term result to hydrostatic / non-thermal-pressure bias and lensing calibration
import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import minimize_scalar,brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    cd={}; exec(open('kernel_scenarios.py').read().split("# ---------- check NFW R500")[0],cd)
    exec(open('stars_test.py').read())
    AL_=dict(z=z.copy(),r25=r25.copy(),fo=fo.copy(),M25=M25.copy(),M5=M5.copy(),xi=xi.copy())
    exec(open('v4_indep.py').read().split("# ---- Vikhlinin 2006")[0])
    GR_=None
G,c,Msun,Mpc=cd['G'],cd['c'],cd['Msun'],cd['Mpc']; lam=cd['lam']; A=cd['A']
R_=lambda a: float(np.sqrt(np.nanmean(np.square(a))))
D=cd['build'](radii_mpc=[.5,1,1.5])
w=D['w']; xi0=D['xi']; Mr0=D['Mr']; M50=D['M5']; fb=D['fb']; g0=D['g']; gd=D['gd']
alpha=.355*(1-np.exp(-xi0/.052))
def xcop(b500,kap):
    """HSE masses are (1-b(w)) x true; b(w)=b500*w^kap (non-thermal fraction grows outward). RGD applied to true masses."""
    bw=np.clip(b500*w**kap,-0.5,0.8); Mt=Mr0/(1-bw); b5=np.clip(b500,-.5,.8); M5t=M50/(1-b5)
    g=g0/(1-bw); xi=xi0*(1/(1-b5)); dens=(Mt/M5t)/w**3; x=fb*g/gd
    K=lambda n,s: 1+lam*np.exp(-x/s)/(-np.expm1(-x/s))**n
    al=.355*(1-np.exp(-xi/.052))
    r_none=np.log10(K(.5+A*xi,1.)*fb); r_v4=np.log10(K((.5+A*xi)*w**(-al),1.)*fb)
    f=lambda p: R_(np.log10(K(.5+A*xi,dens**p)*fb)); o=minimize_scalar(f,bounds=(0,2),method="bounded")
    return R_(r_none),R_(r_v4),o.x,o.fun,f(0.5)
print("A. X-COP 39 pts: bias b(w)=b500*w^kappa applied to the 'true' mass (M_HSE=(1-b)M_true); kappa=0 constant, kappa=1.0 growing outward (non-thermal pressure)")
print("   (b>0: HSE underestimates; Eckert+2019 measure ~6%% at r500 rising outward; b<0 would mean HSE overestimates)")
print("   %6s %5s | rms none  v4   | v5: best p  rms(best)  rms(p=.5)"%("b500","kap"))
for kap in (0.0,1.0):
    for b in (-0.1,0.0,0.06,0.1,0.15,0.2,0.3):
        rn,r4,p,rb,r5=xcop(b,kap); print("   %6.2f %5.1f | %.3f  %.3f | %.2f  %.3f  %.3f"%(b,kap,rn,r4,p,rb,r5))
# ---------------- measured f_gas samples: Vikhlinin (ratio test) and Allen (absolute)
Hv=lambda z_: 71e3/Mpc*np.sqrt(.3*(1+z_)**3+.7)
V=[("A133",.0569,998,1.24,3.14,.065,.084),("A383",.1883,956,1.68,3.10,.090,.122),("A478",.0881,1359,4.23,7.83,.096,.118),("A907",.1603,1117,2.30,4.71,.090,.121),
("A1413",.1429,1339,3.08,7.78,.092,.104),("A1795",.0622,1283,2.66,6.57,.088,.098),("A1991",.0592,753,.60,1.28,.068,.099),("A2029",.0779,1380,4.41,8.29,.090,.121),
("A2390",.2302,1448,3.50,10.88,.127,.138),("MKW4",.0199,635,.30,.74,.044,.063)]
kpc=3.0857e19; Ms=Msun
def pred(M,r,mode,p,P):
    g=G*M/r**2; s=P['dd']**p if mode=="v5" else 1.0; xi_=P['xi']
    nn=(.5+A*xi_)*((r/P['R5'])**(-.355*(1-np.exp(-xi_/.052)))) if mode=="v4" else .5+A*xi_
    try: return brentq(lambda f: f*(1+lam*np.exp(-f*g/P['gd']/s)/(-np.expm1(-f*g/P['gd']/s))**nn)-1,1e-6,5.)
    except ValueError: return np.nan
def vik(b25,b5):
    out=[]
    for n,zz,r5,M25,M5_,f25,f5 in V:
        rc=3*Hv(zz)**2/(8*np.pi*G); R5=r5*kpc; M5h=M5_*1e14*Ms; M25h=M25*1e14*Ms; r25_=(M25h/(2500*rc*4*np.pi/3))**(1/3)
        M5t=M5h/(1-b5); M25t=M25h/(1-b25); gdv=c*Hv(zz)/(2*np.pi); xi_=G*M5t/R5**2/gdv; st=Mstar_tot(M5_*1e14)*Ms
        # gas mass is measured independently of HSE; f_gas = Mgas/M_true
        fb5=(f5*M5h+st)/M5t; fb25=(f25*M25h+.7*st)/M25t
        dd=(M25t/M5t)/(r25_/R5)**3
        P=dict(xi=xi_,gd=gdv,R5=R5,dd=dd)
        out.append((fb25,fb5,pred(M25t,r25_,"none",0,P),pred(M25t,r25_,"v4",0,P),(lambda p,M=M25t,r=r25_,P=P: pred(M,r,"v5",p,P)),pred(M5t,R5,"none",0,P)))
    return out
print("\nB. Vikhlinin: HSE bias applied as (1-b25) at r2500, (1-b500) at r500; f_gas=Mgas/M_true. Ratio test mean ln[(f(2500)/f(500))_obs/(pred)]:")
print("   %5s %5s | none    v4(n)    v5 p=.5  | best p for ratio test"%("b25","b500"))
for b25,b5 in ((0,0),(0.05,0.05),(0.1,0.06),(0.15,0.1),(0.2,0.1),(0.25,0.15),(0.3,0.15),(-0.05,-0.05)):
    o=vik(b25,b5)
    rat=lambda mode,p=.5: np.array([np.log((a[0]/a[1])/((a[2] if mode=="none" else a[3] if mode=="v4" else a[4](p))/a[5])) for a in o])
    f=lambda p: R_(rat("v5",p)); op=minimize_scalar(f,bounds=(0,1.5),method="bounded")
    print("   %5.2f %5.2f | %+.2f   %+.2f(%d)  %+.2f±%.2f | p=%.2f"%(b25,b5,np.mean(rat("none")),np.nanmean(rat("v4")),np.isnan(rat("v4")).sum(),np.mean(rat("v5")),np.std(rat("v5"),ddof=1)/np.sqrt(10),op.x))
# ---------------- Allen: absolute f_gas at r2500, constant b and mass-dependent b
aZ,aR,aF,aM,aM5,aX=AL_['z'],AL_['r25'],AL_['fo'],AL_['M25'],AL_['M5'],AL_['xi']
lMa=np.log10(aM/Msun)
def allen(b0,b1,p,phi=.7):
    b=np.clip(b0+b1*(lMa-14.4),-.5,.8); Mt=aM/(1-b); M5t=aM5/(1-b); xt=aX/(1-b)
    fbo=(aF*aM+phi*Mstar_tot(M5t/Msun)*Msun)/Mt
    def fs(m,r,x_,zz):
        try: return fsolve(m,r,x_,c_*Hz(zz)/(2*np.pi),p)
        except ValueError: return np.nan
    pr=np.array([fs(m,r,x_,zz) for m,r,x_,zz in zip(Mt,aR,xt,aZ)])
    return np.log(fbo/pr)
print("\nC. Allen (42, absolute f_bar at r2500): constant bias b0")
print("   %5s | p=.45: mean  rms | p=.50: mean  rms | best p"%"b0")
for b0 in (-0.1,0,0.1,0.2):
    r45=allen(b0,0,.45); r50=allen(b0,0,.5); op=minimize_scalar(lambda p: R_(allen(b0,0,p)),bounds=(0.1,1.2),method="bounded")
    print("   %5.2f |      %+.3f %.3f |      %+.3f %.3f | %.2f"%(b0,r45.mean(),R_(r45),r50.mean(),R_(r50),op.x))
# ---------------- lensing calibration: CLASH M_L/M_X(2500) = 1.21 median, Mantz 0.96 -> b = 1-M_X/M_true
print("\nD. Lens calibration -> b25 = 1 - 1/(M_L/M_X): Mantz+2016 0.96 -> b=%.2f ; CLASH-overlap median 1.21 -> b=%.2f ; Planck-like 'b=0.2' -> M_L/M_X=%.2f"%(1-1/.96,1-1/1.21,1/.8))
