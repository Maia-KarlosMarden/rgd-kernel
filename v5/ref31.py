# Referee point 3.1: independence, parameter provenance, cross-sample frozen-p prediction with bootstrap uncertainties
import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import brentq,minimize_scalar
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open('stars_test.py').read())      # groups(grows,...), Allen arrays, Mstar_tot, constants
    AL_=dict(z=z.copy(),r25=r25.copy(),fo=fo.copy(),M25=M25.copy(),M5=M5.copy(),xi=xi.copy()); GR_=list(grows)
    exec(open('v4_indep.py').read().split("# ---- Vikhlinin 2006")[0])  # solve(modes none/v4/v5)
    cd={}; exec(open('kernel_scenarios.py').read().split("# ---------- check NFW R500")[0],cd)
rng=np.random.default_rng(11)
PG=np.round(np.linspace(0.0,1.6,33),2)       # p grid
def sol(M,r,xi,gd,mode,w,dd,p): return solve(M,r,xi,gd,mode,w=w,dd=dd,p=p)
# ---------------- datasets: each -> dict(name, n, fb(obs), M,r,xi,gd,w,dens, stars-included)
DS={}
# Vikhlinin 2006 (10 clusters with both radii) at r2500, phi=0.7
Hv=lambda z: 71e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
V=[("A133",.0569,998,1.24,3.14,.065,.084),("A383",.1883,956,1.68,3.10,.090,.122),("A478",.0881,1359,4.23,7.83,.096,.118),("A907",.1603,1117,2.30,4.71,.090,.121),
("A1413",.1429,1339,3.08,7.78,.092,.104),("A1795",.0622,1283,2.66,6.57,.088,.098),("A1991",.0592,753,.60,1.28,.068,.099),("A2029",.0779,1380,4.41,8.29,.090,.121),
("A2390",.2302,1448,3.50,10.88,.127,.138),("MKW4",.0199,635,.30,.74,.044,.063)]
L=[]
for n,zz,r5,M25,M5_,f25,f5 in V:
    rc=3*Hv(zz)**2/(8*np.pi*Gc); R5=r5*kpc; M5k=M5_*1e14*Ms; M25k=M25*1e14*Ms; r25_=(M25k/(2500*rc*4*np.pi/3))**(1/3)
    gd=cc*Hv(zz)/(2*np.pi); xi_=Gc*M5k/R5**2/gd; st=Mstar_tot(M5_*1e14)*Ms
    L.append(dict(M=M25k,r=r25_,xi=xi_,gd=gd,w=r25_/R5,dd=(M25k/M5k)/(r25_/R5)**3,fb=f25+.7*st/M25k,fb500=f5+st/M5k,M5=M5k,R5=R5))
DS["Vikhlinin r2500 (10)"]=L
# Ettori 2019 r200 (9)
H7=lambda zz: 70e3/Mpc*np.sqrt(.3*(1+zz)**3+.7)
X=[("A1644",.0473,3.48,6.69,.128,.126),("A1795",.0622,4.63,6.53,.139,.144),("A2029",.0773,8.65,12.25,.141,.152),("A2142",.0909,8.95,13.64,.158,.168),("A2255",.0809,5.26,10.33,.153,.146),
("A3158",.0597,4.26,6.63,.145,.155),("A644",.0704,5.66,7.67,.132,.139),("A85",.0555,5.65,8.50,.150,.159),("RXC1825",.0650,4.08,6.15,.133,.155)]
L=[]
for n,zz,M5_,M2_,f5,f2 in X:
    rc=3*H7(zz)**2/(8*np.pi*Gc); M5k=M5_*1e14*Ms; M2k=M2_*1e14*Ms
    R5=(M5k/(500*rc*4*np.pi/3))**(1/3); R2=(M2k/(200*rc*4*np.pi/3))**(1/3); gd=cc*H7(zz)/(2*np.pi); xi_=Gc*M5k/R5**2/gd; st=Mstar_tot(M5_*1e14)*Ms
    L.append(dict(M=M2k,r=R2,xi=xi_,gd=gd,w=R2/R5,dd=0.4,fb=f2+st/M2k,fb500=f5+st/M5k,M5=M5k,R5=R5))
DS["X-COP Ettori r200 (9)"]=L
# Groups r2500 (43), phi=0.7
L=[]
for r in GR_:
    st=.7*Mstar_tot(r['M5']/Msun)*Msun
    L.append(dict(M=r['M25'],r=r['R25'],xi=r['xi'],gd=r['gd'],w=r['R25']/r['R5'],dd=5.0*1.0,fb=r['f25']+st/r['M25'],fb500=r['f5']+st/r['M5'],M5=r['M5'],R5=r['R5']))
DS["Groups r2500 (43)"]=L
# Allen 2008 r2500 (42), phi=0.7
L=[]
rc0=lambda zz: 3*Hz(zz)**2/(8*np.pi*G_)
aZ,aR,aF,aM,aM5,aX=AL_['z'],AL_['r25'],AL_['fo'],AL_['M25'],AL_['M5'],AL_['xi']
for i in range(len(aZ)):
    R5=(aM5[i]/(500*rc0(aZ[i])*4*np.pi/3))**(1/3)
    L.append(dict(M=aM[i],r=aR[i],xi=aX[i],gd=c_*Hz(aZ[i])/(2*np.pi),w=aR[i]/R5,dd=5.0,fb=aF[i]+.7*Mstar_tot(aM5[i]/Msun)*Msun/aM[i],fb500=np.nan,M5=aM5[i],R5=R5))
DS["Allen r2500 (42)"]=L
# --------------- residual tables
def res(L,mode,p):
    out=[]
    for d in L:
        pr=solve(d['M'],d['r'],d['xi'],d['gd'],mode,w=d['w'],dd=d['dd'],p=p)
        out.append(np.log(d['fb']/pr) if np.isfinite(pr) else np.nan)
    return np.array(out)
TAB={}
for k,L in DS.items():
    TAB[k]={'none':res(L,'none',0),'v4':res(L,'v4',0)}
    TAB[k]['v5']=np.array([res(L,'v5',p) for p in PG])    # [p, cluster]
# X-COP 39 pts
D0=cd['build'](radii_mpc=[.5,1,1.5]); N=D0['xi'].shape[0]
x=D0['fb']*D0['g']/D0['gd']; xiX=D0['xi']; wX=D0['w']; densX=(D0['Mr']/D0['M5'])/wX**3
Kx=lambda n,s: 1+lam*np.exp(-x/s)/(-np.expm1(-x/s))**n
alpha=.355*(1-np.exp(-xiX/.052))
rX={'none':np.log10(Kx(.5+A*xiX,1.)*D0['fb']),'v4':np.log10(Kx((.5+A*xiX)*wX**(-alpha),1.)*D0['fb'])}
rX['v5']=np.array([np.log10(Kx(.5+A*xiX,densX**p)*D0['fb']) for p in PG])
R_=lambda a: float(np.sqrt(np.nanmean(np.square(a))))
ip=lambda p: int(np.argmin(abs(PG-p)))
print("="*100); print("A. X-COP 39 points (13 clusters x 3 radii): cluster bootstrap (4000) of rms [dex] and of paired differences")
mods={'none (n(xi) at w=1)':lambda ix: R_(rX['none'][ix]),'v4 (2 extra par.)':lambda ix: R_(rX['v4'][ix]),'v5 p=0.50 fixed':lambda ix: R_(rX['v5'][ip(.5)][ix])}
pb=[];bs={k:[] for k in mods}; bsfit=[]
for _ in range(4000):
    ix=rng.integers(0,N,N)
    for k,f in mods.items(): bs[k].append(f(ix))
    v=[R_(rX['v5'][j][ix]) for j in range(len(PG))]; j=int(np.argmin(v)); pb.append(PG[j]); bsfit.append(v[j])
for k,f in mods.items(): a=np.array(bs[k]); print("  %-22s rms %.3f  95%% CI [%.3f, %.3f]"%(k,f(slice(None)),*np.percentile(a,[2.5,97.5])))
a=np.array(bsfit); print("  v5 p fitted            rms %.3f  95%% CI [%.3f, %.3f];  p = %.2f  95%% CI [%.2f, %.2f]"%(R_(rX['v5'][np.argmin([R_(r) for r in rX['v5']])]),*np.percentile(a,[2.5,97.5]),PG[np.argmin([R_(r) for r in rX['v5']])],*np.percentile(pb,[2.5,97.5])))
for a_,b_ in (('v4 (2 extra par.)','v5 p=0.50 fixed'),('none (n(xi) at w=1)','v4 (2 extra par.)')):
    d=np.array(bs[a_])-np.array(bs[b_]); print("  paired rms(%s) - rms(%s): median %+.3f, 95%% CI [%+.3f, %+.3f]; fraction of resamples with v5<v4: %.3f"%(a_,b_,np.median(d),*np.percentile(d,[2.5,97.5]),np.mean(d>0)))
def m2l(r,n): return n*np.log(R_(r)**2)
for Neff,lab in ((39,"N=39 (points independent)"),(13,"N_eff=13 (one independent unit per cluster)")):
    f=Neff/39.
    print("  -2lnL (%s): none %.0f | v4 %.0f (+4) | v5 p=.5 %.0f | v5 fit %.0f (+2);  AIC none/v4/v5(p=.5)/v5(fit) = %.0f / %.0f / %.0f / %.0f"%(lab,
      f*39*np.log(R_(rX['none'])**2),f*39*np.log(R_(rX['v4'])**2),f*39*np.log(R_(rX['v5'][ip(.5)])**2),f*39*np.log(R_(rX['v5'][ip(.45)])**2),
      f*39*np.log(R_(rX['none'])**2),f*39*np.log(R_(rX['v4'])**2)+4,f*39*np.log(R_(rX['v5'][ip(.5)])**2),f*39*np.log(R_(rX['v5'][ip(.45)])**2)+2))
p0=PG[np.argmin([R_(r) for r in rX['v5']])]; print("  X-COP calibrated p (grid) = %.2f -> frozen for everything below"%p0)
# ------------------- B. frozen-p predictions on other samples
print("\n"+"="*100); print("B. Frozen-p (calibrated on X-COP 39 pts only) applied to other samples; no refit. resid = ln(f_bar_obs / f_bar_pred)")
print("   mean = bias, rms = total, sd = scatter about mean; CI from cluster bootstrap (4000)")
for k,L in DS.items():
    n=len(L); print(" %s"%k)
    for lab,r in (("none",TAB[k]['none']),("v4 ",TAB[k]['v4']),("v5 p=%.2f frozen"%p0,TAB[k]['v5'][ip(p0)]),("v5 p=0.50 round",TAB[k]['v5'][ip(.5)])):
        ok=np.isfinite(r); r_=r[ok]; bm=[];br=[]
        for _ in range(4000):
            ix=rng.integers(0,r_.size,r_.size); bm.append(r_[ix].mean()); br.append(np.sqrt(np.mean(r_[ix]**2)))
        print("   %-18s n=%2d  mean %+.3f [%+.3f,%+.3f]  rms %.3f [%.3f,%.3f]  sd %.3f  no-solution %d"%(lab,r_.size,r_.mean(),*np.percentile(bm,[2.5,97.5]),np.sqrt(np.mean(r_**2)),*np.percentile(br,[2.5,97.5]),r_.std(),(~ok).sum()))
# --- r500 control and ratio test
print("\n   Control at r500 (s=1 by construction) and zero-point-free ratio test f_bar(r)/f_bar(r500):")
for k in ("Vikhlinin r2500 (10)","X-COP Ettori r200 (9)"):
    L=DS[k]; c500=np.array([np.log(d['fb500']/solve(d['M5'],d['R5'],d['xi'],d['gd'],'none')) for d in L])
    print("   %s: r500 obs/pred mean %+.3f (SE %.3f)"%(k,c500.mean(),c500.std(ddof=1)/np.sqrt(c500.size)))
    for lab,r in (("none",TAB[k]['none']),("v4 ",TAB[k]['v4']),("v5 .45",TAB[k]['v5'][ip(p0)]),("v5 .50",TAB[k]['v5'][ip(.5)])):
        d=r-c500; ok=np.isfinite(d); d=d[ok]
        print("     ratio test %-4s mean ln[(obs ratio)/(pred ratio)] %+.3f ± %.3f (SE), rms %.3f, n=%d"%(lab,d.mean(),d.std(ddof=1)/np.sqrt(d.size),np.sqrt(np.mean(d**2)),d.size))
# --- C. best-fit p per dataset (heterogeneity)
print("\n"+"="*100); print("C. Best-fit constant p per dataset (95% bootstrap CI): heterogeneity of the radial amplitude")
for k in DS:
    T=TAB[k]['v5']; n=T.shape[1]
    def bestp(ix): v=[R_(T[j][ix]) for j in range(len(PG))]; return PG[int(np.argmin(v))]
    bp=[bestp(rng.integers(0,n,n)) for _ in range(2000)]
    print("   %-26s p = %.2f  [%.2f, %.2f]"%(k,bestp(slice(None)),*np.percentile(bp,[2.5,97.5])))
bp=[]; 
for _ in range(2000):
    ix=rng.integers(0,N,N); v=[R_(rX['v5'][j][ix]) for j in range(len(PG))]; bp.append(PG[int(np.argmin(v))])
print("   %-26s p = %.2f  [%.2f, %.2f]"%("X-COP 39 pts",p0,*np.percentile(bp,[2.5,97.5])))
