# Hierarchical model for the radial amplitude p: p_i ~ N(mu + X_i beta, tau^2 + sigma_i^2); covariates: log M500, z, instrument, log Delta/Delta500
import io,contextlib,numpy as np,warnings,re; warnings.filterwarnings("ignore")
from scipy.optimize import brentq,minimize
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open('stars_test.py').read())
    AL_=dict(z=z.copy(),r25=r25.copy(),fo=fo.copy(),M25=M25.copy(),M5=M5.copy(),xi=xi.copy()); GR_=list(grows); ALL_=list(ALLEN)
Gc=6.674e-11; cc=2.998e8; Ms=1.989e30; Mpc=3.0857e22; kpc=3.0857e19; lam_=.999; A_=.4068
H0=70e3/Mpc; E=lambda z_: np.sqrt(.3*(1+z_)**3+.7)
mf=lambda x: np.log(1+x)-x/(1+x)
def norm(n): return re.sub(r'[^A-Z0-9]','',n.upper().replace('ABELL','A'))
D=[]   # records
# ---- Vikhlinin (Chandra, r2500), phi=0.7
V=[("A133",.0569,998,1.24,3.14,.065,.084),("A383",.1883,956,1.68,3.10,.090,.122),("A478",.0881,1359,4.23,7.83,.096,.118),("A907",.1603,1117,2.30,4.71,.090,.121),
("A1413",.1429,1339,3.08,7.78,.092,.104),("A1795",.0622,1283,2.66,6.57,.088,.098),("A1991",.0592,753,.60,1.28,.068,.099),("A2029",.0779,1380,4.41,8.29,.090,.121),
("A2390",.2302,1448,3.50,10.88,.127,.138),("MKW4",.0199,635,.30,.74,.044,.063)]
Hv=lambda z_: 71e3/Mpc*E(z_)
for n,zz,r5,M25,M5_,f25,f5 in V:
    rc=3*Hv(zz)**2/(8*np.pi*Gc); R5=r5*kpc; M5k=M5_*1e14*Ms; M25k=M25*1e14*Ms; r=(M25k/(2500*rc*4*np.pi/3))**(1/3); st=Mstar_tot(M5_*1e14)*Ms
    D.append(dict(name=n,samp="Vikhlinin",z=zz,instr=0,M=M25k,r=r,M5=M5k,R5=R5,dd=(M25k/M5k)/(r/R5)**3,fb=f25+.7*st/M25k,sf=0.10,H0=71e3/Mpc))
# ---- Ettori 2019 (XMM, r200), all 12 ; stars phi=1
X=[("A1644",.0473,3.48,6.69,.128,.126),("A1795",.0622,4.63,6.53,.139,.144),("A2029",.0773,8.65,12.25,.141,.152),("A2142",.0909,8.95,13.64,.158,.168),("A2255",.0809,5.26,10.33,.153,.146),
("A2319",.0557,7.31,10.18,.189,.237),("A3158",.0597,4.26,6.63,.145,.155),("A3266",.0589,8.80,15.12,.132,.108),("A644",.0704,5.66,7.67,.132,.139),("A85",.0555,5.65,8.50,.150,.159),("RXC1825",.0650,4.08,6.15,.133,.155),("ZW1215",.0766,7.66,13.03,.106,.092)]
for n,zz,M5_,M2_,f5,f2 in X:
    rc=3*H0**2*E(zz)**2/(8*np.pi*Gc); M5k=M5_*1e14*Ms; M2k=M2_*1e14*Ms
    R5=(M5k/(500*rc*4*np.pi/3))**(1/3); R2=(M2k/(200*rc*4*np.pi/3))**(1/3); st=Mstar_tot(M5_*1e14)*Ms
    D.append(dict(name=n,samp="Ettori r200",z=zz,instr=1,M=M2k,r=R2,M5=M5k,R5=R5,dd=0.4,fb=f2+st/M2k,sf=0.10,H0=H0,merger=n in("A2319","A3266","ZW1215")))
# ---- Allen (Chandra, r2500), phi=0.7
rc0=lambda zz: 3*H0**2*E(zz)**2/(8*np.pi*G_)
for i,a in enumerate(ALL_):
    zz=AL_['z'][i]; R5=(AL_['M5'][i]/(500*rc0(zz)*4*np.pi/3))**(1/3)
    D.append(dict(name=a[0],samp="Allen",z=zz,instr=0,M=AL_['M25'][i],r=AL_['r25'][i],M5=AL_['M5'][i],R5=R5,dd=5.0,fb=AL_['fo'][i]+.7*Mstar_tot(AL_['M5'][i]/Msun)*Msun/AL_['M25'][i],sf=a[4]/a[3],H0=H0))
# ---- groups r2500, phi=0.7: Lovisari=XMM, Sun=Chandra
for r in GR_:
    st=.7*Mstar_tot(r['M5']/Msun)*Msun
    D.append(dict(name=r['name'],samp="Groups-L(XMM)" if r['src']=="L" else "Groups-S(Chandra)",z=r['z'],instr=1 if r['src']=="L" else 0,M=r['M25'],r=r['R25'],M5=r['M5'],R5=r['R5'],dd=5.0,fb=r['f25']+st/r['M25'],sf=0.10,H0=r['H']/E(r['z'])))
def fpred(d,p,alpha=1.0,dlnM=0.0):
    M=d['M']*(1+dlnM); gd=cc*d['H0']*E(d['z'])**alpha/(2*np.pi); xi_=Gc*d['M5']*(1+dlnM)/d['R5']**2/gd
    s=d['dd']**p; g=Gc*M/d['r']**2; n=.5+A_*xi_
    h=lambda f: f*(1+lam_*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**n)-1
    try: return brentq(h,1e-6,5.0)
    except ValueError: return np.nan
def solve_p(d,alpha=1.0):
    h=lambda p: np.log(d['fb']/fpred(d,p,alpha))
    try: return brentq(h,-1.0,3.0)
    except Exception: return np.nan
def build(alpha=1.0,sM=0.10,sfmult=1.0):
    rows=[]
    for d in D:
        p=solve_p(d,alpha)
        if not np.isfinite(p): continue
        e=1e-3; dp=(np.log(fpred(d,p+e,alpha))-np.log(fpred(d,p-e,alpha)))/(2*e)
        dM=(np.log(fpred(d,p,alpha,0.01))-np.log(fpred(d,p,alpha,-0.01)))/0.02
        sig=np.sqrt((d['sf']*sfmult)**2+(dM*sM)**2)/abs(dp)
        M5s=d['M5']/Ms
        xi_=Gc*d['M5']/d['R5']**2/(cc*d['H0']*E(d['z'])**alpha/(2*np.pi)); rows.append((d['name'],d['samp'],d['z'],np.log10(M5s)-14.5,d['instr'],np.log10(d['dd']),p,sig,d.get('merger',False),xi_))
    return rows
def ml(rows,cols,fix_tau=None):
    y=np.array([r[6] for r in rows]); s=np.array([r[7] for r in rows]); X=np.column_stack([np.ones(len(rows))]+[np.array([r[c] for r in rows],float) for c in cols])
    def nll(lt):
        t2=np.exp(2*lt) if fix_tau is None else fix_tau**2; w=1/(s**2+t2); A=X.T@(X*w[:,None]); b=np.linalg.solve(A,X.T@(w*y)); r=y-X@b
        return 0.5*np.sum(np.log(s**2+t2)+w*r**2)+0.5*np.log(np.linalg.det(A)) - 0.5*X.shape[1]*0   # REML
    if fix_tau is None:
        from scipy.optimize import minimize_scalar
        o=minimize_scalar(nll,bounds=(-6,1),method="bounded"); lt=o.x; v=o.fun
    else: lt=np.log(fix_tau); v=nll(lt)
    t2=np.exp(2*lt) if fix_tau is None else fix_tau**2; w=1/(s**2+t2); A=X.T@(X*w[:,None]); b=np.linalg.solve(A,X.T@(w*y)); cov=np.linalg.inv(A)
    return dict(tau=np.exp(lt),beta=b,se=np.sqrt(np.diag(cov)),nll=v,n=len(rows),cols=cols)
names={2:"z",3:"logM500-14.5",4:"XMM",5:"log Delta/Delta500"}
def report(tag,rows,cols_list):
    print("\n--- %s  (N=%d)"%(tag,len(rows)))
    for cols in cols_list:
        r=ml(rows,cols); lab=["mu"]+[names[c] for c in cols]
        print("  covariates %-40s tau=%.3f | "%(",".join(names[c] for c in cols) or "none",r['tau'])+"; ".join("%s=%+.3f±%.3f"%(l,b,e) for l,b,e in zip(lab,r['beta'],r['se'])))
if __name__=="__main__":
    rows=build()
    print("N data",len(rows)); 
    for sm in sorted(set(r[1] for r in rows)):
        y=np.array([r[6] for r in rows if r[1]==sm]); s=np.array([r[7] for r in rows if r[1]==sm]); print("%-20s N=%2d  median p_i %.2f  [16,84]%% %.2f..%.2f  median sigma_i %.2f"%(sm,len(y),np.median(y),*np.percentile(y,[16,84]),np.median(s)))
