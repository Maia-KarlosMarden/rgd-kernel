import sys, numpy as np
sys.path.insert(0,'.')
import predict_frozen as pf
from scipy.optimize import brentq
Mpc=pf.Mpc; G=pf.G; Msun=pf.Msun
def Hz(z,H0=70.,Om=0.3): return H0*1e3/Mpc*np.sqrt(Om*(1+z)**3+1-Om)
def rdelta(M,z,Delta=500,H0=70.,Om=0.3):   # Mpc, M in Msun
    rho=3*Hz(z,H0,Om)**2/(8*np.pi*G)
    return (3*M*Msun/(4*np.pi*Delta*rho))**(1/3)/Mpc
def Ez(z): return np.sqrt(0.3*(1+z)**3+0.7)
def fbar(M500,z,r='500',M=None,R=None,**kw):
    R5=rdelta(M500,z)
    if r=='500': return pf.predict_fbar(M500,R5,M500,R5,z,**kw)[0]
    return pf.predict_fbar(M,R,M500,R5,z,**kw)[0]
# ---------- HSC-XXL (2111.10080) Table 1: ln(Y)=alpha+beta*ln(M500 E/1e14)
Mg=(1.95,.08,1.29,(.16,.10)); Ms=(.76,.085,.85,(.12,.09))
def aper(lM,lo,mid,hi):   # quoted corrections at 1e13,1e14,1e15
    return np.interp(lM,[13,14,15],[lo,mid,hi])
def xxl_frac(M,z,rng=None,n=2000):
    Z=np.log(M*Ez(z)/1e14); lM=np.log10(M*Ez(z))
    def one(par,aa,bb):
        a,sa,b,(bu,bd)=par; 
        if rng is None: return np.exp(a+b*Z)*1e12
        da=rng.normal(0,sa,n); db=np.where(rng.random(n)<.5,np.abs(rng.normal(0,bu,n)),-np.abs(rng.normal(0,bd,n)))
        return np.exp(a+da+(b+db)*Z)*1e12
    g=one(Mg,0,0)*aper(lM,1.09,1.04,.97)/Ez(z)/M
    s=one(Ms,0,0)*aper(lM,1.05,1.02,.99)/Ez(z)/M
    return g,s
def run(label,**kw):
    print("\n=== %s  (frozen: %s) ==="%(label,kw))
    rng=np.random.default_rng(1)
    for z in (0.3,):
        print(" HSC-XXL f_b(<r500) vs RGD prediction, z=%.1f"%z)
        print("  M500[1e13]  f_g,obs   f_*,obs   f_b,obs(±)     f_b,RGD   ratio RGD/obs   n_sigma")
        chi=0;N=0
        for M in [1.5e13,3e13,5.6e13,1e14,2e14,3e14,5e14,1e15]:
            g,s=xxl_frac(M,z,rng); fb=g+s; m,sd=fb.mean(),fb.std()
            pr=fbar(M/1e0,z,**kw)
            gm,sm=xxl_frac(M,z)
            print("  %8.1f   %.4f   %.4f   %.4f±%.4f   %.4f    %.3f        %+.2f"%(M/1e13,gm,sm,gm+sm,sd,pr,pr/(gm+sm),(pr-(gm+sm))/sd))
            if 3e13<=M<=5e14: chi+=((pr-(gm+sm))/sd)**2;N+=1
        print("  chi2 (3e13-5e14, %d pts, correlated errs ignored) = %.1f"%(N,chi))
    # z sensitivity at 1e14
    print("  z-sensitivity (M=1e14): "+", ".join("z=%.1f: %.4f"%(z,fbar(1e14,z,**kw)) for z in (.1,.2,.3,.4,.5)))
    # eFEDS binned gas (read by eye from Fig 11) with HSC-XXL stars
    ef=[(3.0,.138,.022,.008),(5.6,.211,.043,.008),(8.4,.263,.052,.006),(13.5,.310,.065,.007),(13.3,.470,.049,.006),(22.5,.697,.053,.008),(25.1,.483,.063,.007),(27.7,.318,.070,.005)]
    print(" eFEDS bins (WL-calibrated M500c; f_gas read by eye ±0.003 systematic): f_gas,RGD = f_bar,RGD - f_*,HSC-XXL")
    print("  bin M500[1e13]  z    f_gas,obs   f_gas,RGD  f_*  resid/sigma")
    tot=0
    for i,(m,z,fg,e) in enumerate(ef,1):
        M=m*1e13; _,s=xxl_frac(M,z); pr=fbar(M,z,**kw); pg=pr-s
        e2=np.hypot(e,.003); tot+=((pg-fg)/e2)**2
        print("  %d  %6.1f  %.2f   %.3f     %.3f     %.3f  %+.2f"%(i,m,z,fg,pg,s,(pg-fg)/e2))
    print("  chi2 = %.1f for 8 points"%tot)
    # Popesso (2411.16555) stacks: eq4 at r500 (nonindependent masses), eq5 at r200
    print(" Popesso+24 relations (masses from Popesso L-M, DM14 NFW; Lovisari M-T: NOT independent of tuning)")
    h=.6774; z=0.15
    def c200(M):
        a=.520+(.905-.520)*np.exp(-.617*z**1.21); b=-.101+.026*z
        return 10**(a+b*np.log10(M/(1e12/h)))
    def nfw_m(x): return np.log(1+x)-x/(1+x)
    print("  M200[1e13]  f_gas500 obs  RGD   | f_gas200 obs  RGD(p=%s)"%kw.get('p',0.5))
    for M200 in [1e13,3e13,1e14,3e14,1e15]:
        R200=rdelta(M200,z,200,67.74,.3089); c=c200(M200)
        r5=brentq(lambda r: nfw_m(c*r/R200)/nfw_m(c)*M200-500/200*(r/R200)**3*M200,0.2*R200,R200)
        # M(<r) = 500/200 * M200*(r/R200)^3 at r500 -> M500
        M500=M200*nfw_m(c*r5/R200)/nfw_m(c)
        R500=r5
        o500=2.23e-7*M500**.39; o200=2.09e-6*M200**.33
        st=pf.Mstar_tot(M500)
        k=dict(kw)
        b5=pf.predict_fbar(M500,R500,M500,R500,z,**k)[0]-st/M500
        b2=pf.predict_fbar(M200,R200,M500,R500,z,**k)[0]-st/M200
        print("  %8.1f    %.4f   %.4f  |  %.4f      %.4f"%(M200/1e13,o500,b5,o200,b2))
if __name__=="__main__":
    run("ARCHIVED prescription (pre-registered)")
    run("Sensitivity: xi_c=0.157, m=8",xic=0.157,m=8)
    run("Sensitivity: p=0.4",p=0.4)
    run("Sensitivity: p=0.6",p=0.6)
