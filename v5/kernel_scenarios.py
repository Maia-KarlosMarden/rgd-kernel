import numpy as np, warnings
from scipy.optimize import least_squares, brentq
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
G=6.674e-11; c=2.998e8; Msun=1.989e30; Mpc=3.0857e22; kpc=3.0857e19
H0=70e3/Mpc; Om=.3; OL=.7; lam=.999; A=.4068; X0=.052; FCOS=0.152
# name,z,c200,M200[1e14],R500,M500
T=[("A85",.0555,3.31,8.50),("A644",.0704,5.58,7.67),("A1644",.0473,1.46,6.69),("A1795",.0622,4.55,6.53),
("A2029",.0773,4.26,12.57),("A2142",.0909,3.14,13.64),("A2255",.0809,1.37,10.33),("A2319",.0557,4.86,10.18),
("A3158",.0597,2.88,6.63),("A3266",.0589,2.04,15.12),("HydraA",.0538,5.51,3.01),("RXC1825",.0650,3.35,6.15),("ZW1215",.0766,2.11,13.03)]
TAB_R500=np.array([1.235,1.230,1.054,1.153,1.423,1.424,1.196,1.346,1.123,1.430,0.904,1.105,1.358])
N=len(T); mf=lambda x: np.log(1+x)-x/(1+x)
Hz=lambda z: H0*np.sqrt(Om*(1+z)**3+OL)
rms=lambda a: float(np.sqrt(np.mean(np.asarray(a)**2)))
def Kern(kind,p,x,xi,w):
    nxi=0.5+A*xi; sat=1-np.exp(-xi/X0); amp=1.0; s=1.0; n=nxi
    if kind=="M1": n=nxi*w**(-p[0]*sat)
    elif kind=="M2": s=w**(-p[0]*sat)
    elif kind=="M3": amp=w**(-p[0]*sat)
    elif kind=="M6": n=nxi*w**(-p[0]*sat); s=w**(-p[1]*sat)
    xs=x/s; q=np.exp(-xs); om=-np.expm1(-xs)
    return 1+lam*amp*q/om**n
def build(b=0.0,kap=None,fshift=0.0,fmult=1.0,radii_mpc=None,radii_w=None):
    D=dict(g=[],w=[],xi=[],fb=[],gd=[],M=[],c=[],z=[],M5=[],R5=[],Mr=[]);
    for (nm,z,cc,M2) in T:
        rho=3*Hz(z)**2/(8*np.pi*G); gd=c*Hz(z)/(2*np.pi)
        M200=M2*1e14*Msun; R200=(3*M200/(800*np.pi*rho))**(1/3)
        Mf=lambda r: M200*mf(cc*r/R200)/mf(cc)/(1-b)
        R5=brentq(lambda R: Mf(R)-500*4*np.pi/3*rho*R**3,0.2*Mpc,3*Mpc); M5=Mf(R5); xi=G*M5/R5**2/gd
        r=np.array(radii_mpc)*Mpc if radii_mpc is not None else np.array(radii_w)*R5
        w=r/R5; Mr=Mf(r); g=G*Mr/r**2
        fg5=0.0708*(M5/Msun/1e13)**0.135
        if kap is None: fb=np.full_like(w,fg5*1.09)
        else:
            Mg5=fg5*M5; Ms5=0.39*(M5/Msun/1e14)**-0.84*Mg5
            fb=fg5*w**kap+Ms5*mf(2.9*w)/mf(2.9)/Mr
        fb=fb*fmult+fshift
        for k,v in dict(g=g,w=w,xi=np.full_like(w,xi),fb=fb,gd=np.full_like(w,gd),Mr=Mr,c=np.full_like(w,cc),z=np.full_like(w,z),M5=np.full_like(w,M5),R5=np.full_like(w,R5)).items(): D[k].append(v)
    return {k:np.array(v) for k,v in D.items() if k!="M"}
def resid(kind,p,D):
    x=D['fb']*D['g']/D['gd']; return np.log10(Kern(kind,p,x,D['xi'],D['w'])*D['fb'])
P0={"M1":[.35],"M2":[.8],"M3":[.6],"M6":[.35,.5]}
def fit(kind,D,idx=None):
    sl=slice(None) if idx is None else idx
    f=least_squares(lambda p: resid(kind,p,D)[sl].ravel(),P0[kind],bounds=(-2,5)); return f.x
# ---------- check NFW R500 vs table
D0=build(radii_mpc=[.5,1,1.5]); print("check R500(NFW) vs table: max rel diff = %.3f%%"%(100*np.max(np.abs(D0['R5'][:,0]/Mpc/TAB_R500-1))))
KINDS=["M1","M2","M3"]
print("\n=== S0 baseline (paper mode: f_bar=1.09 f_gas const; 39 pts) ===")
for k in KINDS+["M6"]:
    p=fit(k,D0); r=resid(k,p,D0); print(" %s p=%s RMS=%.4f per-radius %s"%(k,np.round(p,3),rms(r),[round(rms(r[:,j]),3) for j in range(3)]))
print(" M0 (n(xi) only, w=1 law): RMS=%.4f"%rms(np.log10(Kern("M0",[0],D0['fb']*D0['g']/D0['gd'],D0['xi'],D0['w'])*D0['fb'])))

print("\n=== S5 cluster bootstrap (4000 resamples, paper mode) ===")
grid=np.linspace(-.5,2.5,301); R={k:np.array([resid(k,[pp],D0) for pp in grid]) for k in ("M1","M2","M3")}
rng=np.random.default_rng(1); res={k:[] for k in R}; pbest={k:[] for k in R}
for _ in range(4000):
    ix=rng.integers(0,N,N)
    for k in R:
        v=np.sqrt(np.mean(R[k][:,ix,:]**2,axis=(1,2))); j=v.argmin(); res[k].append(v[j]); pbest[k].append(grid[j])
res={k:np.array(v) for k,v in res.items()}; pbest={k:np.array(v) for k,v in pbest.items()}
for k in R: print(" %s best-param %.2f [%.2f,%.2f] 95%%; RMS %.4f [%.4f,%.4f]"%(k,np.median(pbest[k]),*np.percentile(pbest[k],[2.5,97.5]),np.median(res[k]),*np.percentile(res[k],[2.5,97.5])))
for a,b2 in (("M1","M2"),("M1","M3"),("M2","M3")):
    d=res[a]-res[b2]; print(" P(%s has lower RMS than %s) = %.3f ; dRMS(%s-%s)=%.4f [%.4f,%.4f]"%(b2,a,np.mean(d>0),a,b2,np.median(d),*np.percentile(d,[2.5,97.5])))

print("\n=== S6 per-cluster fit of the radial parameter & correlations (paper mode) ===")
pg=np.linspace(-.5,3,351); out={}
for k in ("M1","M2"):
    Rk=np.array([resid(k,[pp],D0) for pp in pg]); pc=np.array([pg[np.argmin(np.sqrt(np.mean(Rk[:,i,:]**2,axis=1)))] for i in range(N)]); out[k]=pc
    cc_=D0['c'][:,0]; xi_=D0['xi'][:,0]; M5=D0['M5'][:,0]/Msun
    print(" %s per-cluster param: median %.2f, range [%.2f,%.2f]; Spearman vs c200 %.2f, xi %.2f, M500 %.2f"%(k,np.median(pc),pc.min(),pc.max(),spearmanr(pc,cc_)[0],spearmanr(pc,xi_)[0],spearmanr(pc,M5)[0]))
for k in ("M1","M2"):
    p=fit(k,D0); r=resid(k,p,D0)[:,0]; print(" %s residual at 0.5 Mpc vs c200: Spearman %.2f (p=%.3f); vs xi %.2f"%(k,spearmanr(r,D0['c'][:,0])[0],spearmanr(r,D0['c'][:,0])[1],spearmanr(r,D0['xi'][:,0])[0]))

print("\n=== S3/S4 robustness: hydrostatic bias, f_bar normalisation (paper mode, refit 1 param each) ===")
print(" scenario                      |  M1: a    RMS  |  M2: b    RMS  |  M3: g    RMS")
for lab,kw in (("baseline",{}),("HSE bias b=0.10",dict(b=.10)),("HSE bias b=0.20",dict(b=.20)),("f_bar -0.02",dict(fshift=-.02)),("f_bar +0.02",dict(fshift=.02)),("f_bar x0.85",dict(fmult=.85)),("f_bar x1.15",dict(fmult=1.15))):
    D=build(radii_mpc=[.5,1,1.5],**kw); s=" %-29s |"%lab
    for k in KINDS:
        p=fit(k,D); s+="  %5.2f %.4f |"%(p[0],rms(resid(k,p,D)))
    print(s)

print("\n=== S1 inner-radius extrapolation (profile mode f_bar=f_gas(w^kap)+f_*(NFW c=2.9)); calibrate on 0.5/1/1.5 Mpc, predict w=0.1-0.3 ===")
Wi=[0.1,0.15,0.2,0.3]; Mods=["M0","M1","M2","M3"]
for kap in (0.15,0.30,0.45):
    Dout=build(kap=kap,radii_mpc=[.5,1,1.5]); Din=build(kap=kap,radii_w=Wi)
    print(" kappa_gas=%.2f  f_bar(R500)~%.3f, f_bar(0.1R500)~%.3f"%(kap,np.mean(build(kap=kap,radii_w=[1.0])['fb']),np.mean(Din['fb'][:,0])))
    for k in Mods:
        p=[0.0] if k=="M0" else fit(k,Dout)
        ro=resid(k,p,Dout) if k!="M0" else np.log10(Kern("M0",p,Dout['fb']*Dout['g']/Dout['gd'],Dout['xi'],Dout['w'])*Dout['fb'])
        ri=resid(k,p,Din) if k!="M0" else np.log10(Kern("M0",p,Din['fb']*Din['g']/Din['gd'],Din['xi'],Din['w'])*Din['fb'])
        print("   %s p=%s in-sample RMS=%.4f | inner bins (w=%s): mean %s  RMS %s"%(k,np.round(p,2),rms(ro),Wi,np.round(ri.mean(axis=0),3),np.round([rms(ri[:,j]) for j in range(len(Wi))],3)))

print("\n=== S2 implied baryon profile from kernel + total mass (profile mode kap=0.30; calibrated on 0.5/1/1.5 Mpc) ===")
kap=.30; Dout=build(kap=kap,radii_mpc=[.5,1,1.5]); Wg=np.geomspace(0.03,2.0,36); Dg=build(kap=kap,radii_w=Wg)
fs=np.geomspace(1e-3,0.95,500)
for k in ("M1","M2","M3"):
    p=fit(k,Dout); nroot0=0; over=0; noroot=0; multi=0; decr=0; tot=0; fm={0.1:[],0.2:[],0.3:[],0.5:[],1.0:[]}
    for i in range(N):
        fimp=[]; 
        for j,wj in enumerate(Wg):
            g=Dg['g'][i,j]; gd=Dg['gd'][i,j]; xi=Dg['xi'][i,j]; Kf=fs*Kern(k,p,fs*g/gd,xi,wj)-1
            idx=np.where(np.diff(np.sign(Kf))!=0)[0]; roots=[brentq(lambda f: f*Kern(k,p,f*g/gd,xi,wj)-1,fs[a],fs[a+1]) for a in idx]
            tot+=1
            if len(roots)==0: noroot+=1; fimp.append(np.nan); continue
            if len(roots)>1: multi+=1
            f=min(roots,key=lambda r_: abs(r_-Dg['fb'][i,j])); fimp.append(f); over+=f>FCOS
        fimp=np.array(fimp); Mb=fimp*Dg['Mr'][i]; ok=~np.isnan(Mb); mb=Mb[ok]; decr+=int(np.sum(np.diff(mb)<-1e-9*mb[:-1]))
        for wv in fm:
            jj=np.argmin(abs(Wg-wv)); fm[wv].append(fimp[jj])
    print(" %s p=%.2f | no solution %d/%d, >1 solution %d/%d, f_implied>f_cosmic(%.3f): %d, decreasing-M_bar steps: %d"%(k,p[0],noroot,tot,multi,tot,FCOS,over,decr))
    print("     median implied f_bar at w=%s : %s   (assumed profile: %s)"%(list(fm),[round(float(np.nanmedian(fm[w_])),3) for w_ in fm],[round(float(np.median(build(kap=kap,radii_w=[w_])['fb'])),3) for w_ in fm]))

print("\n=== S10 XLSSC 122 core (z=1.98, M500=1.66e14; NFW 2.6e14 c=6.3) with gas+stars profile ===")
zc=1.98; rho=3*Hz(zc)**2/(8*np.pi*G); gd=c*Hz(zc)/(2*np.pi); M5=1.66e14*Msun; R5=(3*M5/(4*np.pi*500*rho))**(1/3); xi=G*M5/R5**2/gd
M200=2.6e14*Msun; R200=(3*M200/(800*np.pi*rho))**(1/3); Mc=M200*mf(6.3*100*kpc/R200)/mf(6.3); gc=G*Mc/(100*kpc)**2; wc=100*kpc/R5
fg5=0.0708*(1.66e14/1e13)**0.135; Ms5r=0.39*(1.66)**-0.84
for kap in (0.15,0.30,0.45):
    Dout=build(kap=kap,radii_mpc=[.5,1,1.5]); Mg5=fg5*M5; fb_c=fg5*wc**kap+Ms5r*Mg5*mf(2.9*wc)/mf(2.9)/Mc
    s=" kappa=%.2f f_bar(core)=%.3f :"%(kap,fb_c)
    for k in ("M1","M2","M3"):
        p=fit(k,Dout); Kc=Kern(k,p,fb_c*gc/gd,xi,wc); s+="  %s dlog=%+.3f (K=%.1f)"%(k,np.log10(Kc*fb_c),Kc)
    print(s)
