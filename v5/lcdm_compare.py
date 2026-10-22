import io,contextlib,numpy as np
from scipy.optimize import least_squares, brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open('sparc_mbar.py').read().split("base=np.sqrt")[0])   # SPARC arrays: gb, gobs, ids, R, m, gl, ...
    Gs=G
    exec(open('groups_test.py').read())
mf=lambda x: np.log(1+x)-x/(1+x)
rho_c0=3*(70e3/3.0857e22)**2/(8*np.pi*6.674e-11)
# ---------- SPARC: NFW halos on top of the same baryons (same Upsilon) ----------
def gnfw(r_m,M200,c,rhoc):
    R200=(3*M200/(4*np.pi*200*rhoc))**(1/3); return 6.674e-11*M200*mf(c*r_m/R200)/mf(c)/r_m**2
mm=m&np.isfinite(gb)&(gb>0)
def fit_gal(g,free_c=True):
    i=np.where((ids==g)&mm)[0]; r_m=R[i]*kpc
    if free_c: f=lambda q: np.log10(gb[i]+gnfw(r_m,10**q[0]*Msun,10**q[1],rho_c0))-np.log10(gobs[i]); best=None
    else:
        def f(q):
            M=10**q[0]*Msun; c=5.71*(M/Msun/(2e12/.7))**-0.084; return np.log10(gb[i]+gnfw(r_m,M,c,rho_c0))-np.log10(gobs[i])
    best=None
    for l0 in (9.5,10.5,11.5,12.5):
        for c0 in ((0.6,1.0,1.4) if free_c else (0,)):
            try:
                s=least_squares(f,[l0,c0] if free_c else [l0],bounds=([7,0.3] if free_c else [7],[14,1.8] if free_c else [14]))
            except Exception: continue
            if best is None or s.cost<best.cost: best=s
    return best.fun
gal_ok=[g for g in gl if ((ids==g)&mm).sum()>=5]
res2=np.concatenate([fit_gal(g,True) for g in gal_ok]); res1=np.concatenate([fit_gal(g,False) for g in gal_ok])
sel=np.isin(ids,gal_ok)&mm
rgd=np.log10(gb*K(gb/gd))[sel]-np.log10(gobs[sel])
n=sel.sum(); rm_=lambda a: np.sqrt(np.mean(a**2)); aic=lambda r,k: n*np.log(np.mean(r**2))+2*k
print("SPARC: %d galaxies, %d points (Upsilon_disk=0.5, bulge 0.7 for all models)"%(len(gal_ok),n))
print(" RGD kernel (s=1, g_dagger=cH0/2pi)      : free params 0                rms %.3f dex  AIC %8.0f"%(rm_(rgd),aic(rgd,0)))
print(" LCDM NFW + Duffy c(M)                    : 1 free/galaxy (M200)    = %3d rms %.3f dex  AIC %8.0f"%(len(gal_ok),rm_(res1),aic(res1,len(gal_ok))))
print(" LCDM NFW free (M200,c)                   : 2 free/galaxy           = %3d rms %.3f dex  AIC %8.0f"%(2*len(gal_ok),rm_(res2),aic(res2,2*len(gal_ok))))
# ---------- Groups: predict M(<r2500) ----------
def nfw_from_M500(M5,z,H):
    rc_=3*H**2/(8*np.pi*6.674e-11)
    def M200_of(M200):
        c=5.71*(M200/Msun/(2e12/.7))**-0.084*(1+z)**-0.47; R200=(3*M200/(4*np.pi*200*rc_))**(1/3); R5=R200*0 ; return c,R200
    def eq(lM200):
        M200=10**lM200; c,R200=M200_of(M200)
        # find r500: M(<r)=500 rho_c 4pi/3 r^3
        f=lambda lr: (M200*mf(c*10**lr/R200)/mf(c))-500*rc_*4*np.pi/3*(10**lr)**3
        r5=10**brentq(f,np.log10(R200*0.05),np.log10(R200*1.2)); return M200*mf(c*r5/R200)/mf(c)-M5
    l=brentq(eq,np.log10(M5)-0.5,np.log10(M5)+1.2); M200=10**l; c,R200=M200_of(M200); return M200,c,R200,rc_
resN=[];resR=[]
for r in grows:
    M200,c,R200,rc_=nfw_from_M500(r['M5'],r['z'],r['H'])
    # NFW prediction of M(<r2500): use observed r2500 (defined as mean density 2500 rho_c, so solve NFW own r2500 -> predicts the radius; compare mass at the OBSERVED radius)
    Mp=M200*mf(c*r['R25']/R200)/mf(c); resN.append(np.log10(Mp/r['M25']))
    sat=1.0; stars=0.09; fb=r['f25']*(1+stars); x=fb*6.674e-11*r['M25']/r['R25']**2/r['gd']; s=5.0**(.5*sat); n_=.5+A*r['xi']; xs=x/s
    Kp=1+lam*np.exp(-xs)/(-np.expm1(-xs))**n_; resR.append(np.log10(Kp*fb))
resN=np.array(resN); resR=np.array(resR)
print("\nGROUPS (43): predict mass inside observed r2500. LCDM: NFW+Duffy c(M200) normalised to M500, 0 free params. RGD: switch ON, p=1/2, baryons within r2500 + stars 0.09 (M500 only enters via xi)")
print(" LCDM NFW+Duffy: rms %.3f mean %+.3f | RGD (p=0.5): rms %.3f mean %+.3f | RGD with p=0.8: "%(rm_(resN),resN.mean(),rm_(resR),resR.mean()),end="")
r8=[]
for r in grows:
    fb=r['f25']*1.09; x=fb*6.674e-11*r['M25']/r['R25']**2/r['gd']; s=5.0**(.8); n_=.5+A*r['xi']; xs=x/s; r8.append(np.log10((1+lam*np.exp(-xs)/(-np.expm1(-xs))**n_)*fb))
r8=np.array(r8); print("rms %.3f mean %+.3f"%(rm_(r8),r8.mean()))
