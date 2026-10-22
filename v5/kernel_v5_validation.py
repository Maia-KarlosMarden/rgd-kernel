import numpy as np, warnings, time, sys
exec(open('kernel_scenarios.py').read().split("# ---------- check NFW R500")[0])
from scipy.stats import spearmanr
T0=time.time()
def KG(x,xi,s=1.0,amp=1.0):
    n=0.5+A*xi; xs=x/s; return 1+lam*amp*np.exp(-xs)/(-np.expm1(-xs))**n
def dens(D): return (D['g']/(G*D['M5']/D['R5']**2))/D['w']
def Kmodel(kind,p,x,xi,w,d):
    sat=1-np.exp(-xi/X0)
    if kind in("P","Pfix"): return Kern("M1",[p[0]],x,xi,w)
    if kind=="M2": return Kern("M2",[p[0]],x,xi,w)
    if kind=="M0": return Kern("M0",[0.0],x,xi,w)
    return KG(x,xi,s=d**(p[0]*sat))            # V5 (p fixed 0.5) or V5f (p fitted)
def rm(kind,p,D): return np.log10(Kmodel(kind,p,D['fb']*D['g']/D['gd'],D['xi'],D['w'],dens(D))*D['fb'])
START={"P":0.35,"M2":0.8,"V5f":0.45}
def fitm(kind,D,rows=slice(None),cols=slice(None)):
    if kind=="V5": return [0.5]
    if kind=="Pfix": return [0.355]
    if kind=="M0": return [0.0]
    f=least_squares(lambda p: rm(kind,p,D)[rows][:,cols].ravel(),[START[kind]],bounds=(-2,5)); return list(f.x)
LAB={"Pfix":"paper(a=.355)","P":"paper(a fit)","M2":"cut w^-b","V5":"V5 p=1/2","V5f":"V5 p fit"}
ORDER=["Pfix","P","M2","V5","V5f"]
Dpap=build(radii_mpc=[.5,1,1.5]); SC={}
print("##### A. REGRESSION AT w=1 (calibration radius): all kernels must coincide #####")
D1=build(radii_w=[1.0]); r={k:rm(k,fitm(k,Dpap),D1) for k in ORDER}
print(" X-COP at R500, paper-mode f_bar: RMS=%.4f mean=%+.4f (paper reports 0.015 / +0.006) ; max |difference between kernels| = %.2e dex"%(rms(r["V5"]),r["V5"].mean(),max(np.max(np.abs(r[k]-r["V5"])) for k in ORDER)))
print("\n##### B. STRUCTURAL PROPERTIES of V5 kernel #####")
bad=0; tot=0
for xi in (0.02,0.1,0.3,0.55,0.8):
    for s in (0.05,0.3,1,3,10,60):
        x=np.logspace(-4,1.7,3000); fx=x*KG(x,xi,s=s); tot+=1; bad+=int(not np.all(np.diff(fx)>0)); 
        assert np.all(KG(x,xi,s=s)>=1)
print(" g_pred=x*K(x/s) strictly increasing in g_bar: %d/%d parameter combos fail ; K>=1 always ; max deep-limit exponent n(xi<=0.8)=%.3f (<1)"%(bad,tot,0.5+A*0.8))
print("\n##### C. PHANTOM-MASS POSITIVITY along radius (profile mode kappa=0.30, calibrated on 0.5/1/1.5 Mpc) #####")
kap=.3; Do=build(kap=kap,radii_mpc=[.5,1,1.5]); Wg=np.geomspace(0.03,2.5,60); Dg=build(kap=kap,radii_w=Wg)
for k in ORDER:
    p=fitm(k,Do); Kg=Kmodel(k,p,Dg['fb']*Dg['g']/Dg['gd'],Dg['xi'],Dg['w'],dens(Dg)); Mph=(Kg-1)*Dg['fb']*Dg['Mr']; Mpred=Kg*Dg['fb']*Dg['Mr']
    viol=int(np.sum(np.diff(Mph,axis=1)< -1e-9*Mph[:,:-1])); rr=np.log10(Kg*Dg['fb'])
    print(" %-14s p=%s  decreasing phantom-mass steps: %3d/%d | residual vs NFW at w=0.05/0.1/0.3/1: %s"%(LAB[k],np.round(p,2),viol,Mph.shape[0]*(len(Wg)-1),[round(float(np.median(rr[:,np.argmin(abs(Wg-w_))])),3) for w_ in (0.05,0.1,0.3,1.0)]))
print("\n##### D. INFORMATION CRITERIA on the 39 X-COP points (paper mode) #####")
def nfw_duffy():
    rs=[]
    for i,(nm,z,cc,M2) in enumerate(T):
        rho=3*Hz(z)**2/(8*np.pi*G)
        def res(lm):
            lm=float(np.atleast_1d(lm)[0]); M200=10**lm*Msun; R200=(3*M200/(800*np.pi*rho))**(1/3); ccc=5.71*(M200/Msun/(2e12/.7))**-0.084*(1+z)**-0.47
            r_=np.array([.5,1,1.5])*Mpc; Mtab=Dpap['Mr'][i]; return np.log10(M200*mf(ccc*r_/R200)/mf(ccc)/Mtab)
        f=least_squares(res,[14.5]); rs+=list(res(f.x))
    return np.array(rs)
rn=nfw_duffy(); n_=39; AIC=lambda r_,k: n_*np.log(np.mean(r_**2))+2*k
print(" NFW+Duffy c(M) (13 free M200): RMS=%.4f AIC=%.1f"%(rms(rn),AIC(rn,13)))
for k in ORDER:
    rr=rm(k,fitm(k,Dpap),Dpap); kfit={"Pfix":0,"P":1,"M2":1,"V5":0,"V5f":1}[k]
    print(" %-14s RMS=%.4f | AIC(k=%d fitted)=%.1f ; AIC(+13 inputs M500)=%.1f ; vs NFW dAIC(13 vs 13)=%+.1f"%(LAB[k],rms(rr),kfit,AIC(rr,kfit),AIC(rr,kfit+13),AIC(rr,kfit+13)-AIC(rn,13)))
print("\n##### E. ROBUSTNESS FACTORIAL: kappa x HSE-bias x f_bar shift (36 scenarios, params refit each) #####")
res={k:[] for k in ORDER}; inn={k:[] for k in ORDER}; scen=[]
for kap in (None,0.15,0.30,0.45):
    for b in (0.0,0.1,0.2):
        for fs in (-0.02,0.0,0.02):
            Do=build(b=b,kap=kap,fshift=fs,radii_mpc=[.5,1,1.5]); scen.append((kap,b,fs))
            Di=build(b=b,kap=kap,fshift=fs,radii_w=[0.1,0.2]) if kap is not None else None
            for k in ORDER:
                p=fitm(k,Do); res[k].append(rms(rm(k,p,Do)))
                if Di is not None: inn[k].append(float(np.mean(rm(k,p,Di)[:,0])))
R=np.array([res[k] for k in ORDER]); 
for a,k in enumerate(ORDER): print(" %-14s in-sample RMS over 36 scen: min %.4f median %.4f max %.4f ; inner(w=0.1) mean resid over 27 profile scen: median %+.3f range [%+.3f,%+.3f]"%(LAB[k],R[a].min(),np.median(R[a]),R[a].max(),np.median(inn[k]),min(inn[k]),max(inn[k])))
iV=ORDER.index("V5"); print(" V5(p=1/2,no fit) has lower RMS than: paper(a=.355) in %d/36, paper(a fit) %d/36, cut-w %d/36 scenarios"%tuple(int(np.sum(R[iV]<R[ORDER.index(k)])) for k in ("Pfix","P","M2")))
print("\n##### F. WITHIN-DATA RADIAL VALIDATION (calibrate on a subset of radii, predict the held-out radius) #####")
for lab,kw in (("paper f const",{}),("profile k=.15",dict(kap=.15)),("profile k=.30",dict(kap=.30)),("profile k=.45",dict(kap=.45))):
    Dm=build(radii_mpc=[.5,1,1.5],**kw); s_=" %-14s predict 0.5Mpc (fit on 1.0+1.5):"%lab
    for k in ORDER:
        p=fitm(k,Dm,cols=[1,2]); r_=rm(k,p,Dm)[:,0]; s_+=" %s %.3f/%+.3f;"%({"Pfix":"Pf","P":"P","M2":"M2","V5":"V5","V5f":"V5f"}[k],rms(r_),r_.mean())
    print(s_)
for lab,kw in (("paper f const",{}),("profile k=.30",dict(kap=.30))):
    Dm=build(radii_mpc=[.5,1,1.5],**kw); s_=" %-14s predict 1.5Mpc (fit on 0.5+1.0):"%lab
    for k in ORDER:
        p=fitm(k,Dm,cols=[0,1]); r_=rm(k,p,Dm)[:,2]; s_+=" %s %.3f/%+.3f;"%({"Pfix":"Pf","P":"P","M2":"M2","V5":"V5","V5f":"V5f"}[k],rms(r_),r_.mean())
    print(s_)
print("   (format: RMS dex / mean residual dex)")
print("\n##### G. CLUSTER BOOTSTRAP (4000): V5 (p=1/2, nothing fitted) vs kernels refit on each resample #####")
grid=np.linspace(-.5,2.5,301); RG={k:np.array([rm(k,[q],Dpap) for q in grid]) for k in ("P","M2")}; rV=rm("V5",[0.5],Dpap)
rng=np.random.default_rng(3); dP=[];dM=[]
for _ in range(4000):
    ix=rng.integers(0,N,N); v=rms(rV[ix]); dP.append(np.sqrt(np.mean(RG["P"][:,ix,:]**2,axis=(1,2))).min()-v); dM.append(np.sqrt(np.mean(RG["M2"][:,ix,:]**2,axis=(1,2))).min()-v)
dP=np.array(dP); dM=np.array(dM); print(" P(V5 lower RMS than paper-refit)=%.3f dRMS %.4f [%.4f,%.4f] ; P(V5 lower than cut-w refit)=%.3f dRMS %.4f [%.4f,%.4f]"%(np.mean(dP>0),np.median(dP),*np.percentile(dP,[2.5,97.5]),np.mean(dM>0),np.median(dM),*np.percentile(dM,[2.5,97.5])))
print("\n##### H. WHAT STRUCTURE IS LEFT IN V5 RESIDUALS (paper mode) #####")
rV=rm("V5",[0.5],Dpap); cc=Dpap['c'][:,0]
for nm_,v in (("c200",cc),("xi",Dpap['xi'][:,0]),("M500",Dpap['M5'][:,0]),("z",Dpap['z'][:,0])): print(" Spearman(resid@0.5Mpc, %s) = %+.2f ; Spearman(cluster-mean resid, %s) = %+.2f"%(nm_,spearmanr(rV[:,0],v)[0],nm_,spearmanr(rV.mean(1),v)[0]))
def r_pc(p,D):
    sat=1-np.exp(-D['xi']/X0); s=dens(D)**((p[0]+p[1]*(D['c']-3.5))*sat); return np.log10(KG(D['fb']*D['g']/D['gd'],D['xi'],s=s)*D['fb'])
lo1=[];lo2=[]
for i in range(N):
    keep=[j for j in range(N) if j!=i]
    q1=least_squares(lambda p: rm("V5f",p,Dpap)[keep].ravel(),[.45]).x; q2=least_squares(lambda p: r_pc(p,Dpap)[keep].ravel(),[.45,0.0]).x
    lo1+=list(rm("V5f",q1,Dpap)[i]); lo2+=list(r_pc(q2,Dpap)[i])
q2=least_squares(lambda p: r_pc(p,Dpap).ravel(),[.45,0.]).x
print(" LOCO RMS: V5f(1 param)=%.4f ; V5 with p=a+b(c200-3.5) (2 params)=%.4f  [a=%.3f b=%+.3f] -> residual concentration signal ~ %.0f%% of V5 error"%(rms(lo1),rms(lo2),q2[0],q2[1],100*(1-rms(lo2)/rms(lo1))))
np.save('Dpap_dummy.npy',np.array([0]))
print("\nT=%.0fs"%(time.time()-T0))
