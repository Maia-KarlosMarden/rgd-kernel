"""PRE-REGISTERED (written before running). Search for an equation for p.
Data: per-system p_i (solve f_bar(p)=f_obs) for Vikhlinin, Ettori r200, Allen, groups (ref41 machinery) + XXL-365 C1 (162 systems, r=500 kpc).
Hierarchical ML: p_i ~ N(X_i b, tau^2 + sigma_i^2).  Candidate covariates (available for every system): z, logM500, log Delta/Delta500 at the measurement radius (ldd),
  log10 xi (= g500/g_dagger), XMM-vs-Chandra dummy.  Models: constant; each single covariate; (logM,z); (logM,ldd); (logxi,ldd); (logM,z,ldd); + instrument dummy.
Selection: ML (not REML) AIC = 2 nll + 2 (k+1), AND leave-one-sample-out (LOSO) predictive log score for the held-out sample (sum over its systems of log N(p_i; X b, tau^2+s_i^2)).
Extra XXL-only exploration (declared): T300kpc, Ngal, Lx residual vs M-T mass, r_ap/r500.
A model is only reported as an 'equation' if it beats constant p in LOSO score summed over samples AND in AIC.  Stars: Kravtsov (archived).
"""
import io,contextlib,sys,numpy as np,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import ref41 as R
    import blind_xxl365_test as X
D=R.D; Ms=R.Ms; Mpc=R.Mpc
# ---- add XXL C1
extra={}
for d in X.out:
    if d['cls']!=1: continue
    Mr=d['Mr']*Ms; M5=d['M5']*Ms; R5=d['r5']*Mpc; r=0.5*Mpc
    st=R.Mstar_tot(d['M5'])*Ms
    fb=d['Mg']*Ms/Mr+st/Mr; sf=np.hypot(d['eMg']*1e11/d['Mg'],0.15)
    nm="XLSSC%d"%d['id']
    D.append(dict(name=nm,samp="XXL C1",z=d['z'],instr=1,M=Mr,r=r,M5=M5,R5=R5,dd=(Mr/M5)/(r/R5)**3,fb=fb,sf=sf,H0=R.H0))
    extra[nm]=d
rows=R.build()
rows=[r for r in rows if np.isfinite(r[6]) and r[7]>0 and -1<r[6]<3]
print("N systems with p_i:",len(rows)," by sample:",{s:sum(1 for r in rows if r[1]==s) for s in sorted(set(r[1] for r in rows))})
y=np.array([r[6] for r in rows]); s=np.array([r[7] for r in rows]); samp=np.array([r[1] for r in rows])
cov=dict(z=np.array([r[2] for r in rows]),logM=np.array([r[3] for r in rows]),ldd=np.array([r[5] for r in rows]),
         lxi=np.log10(np.array([r[9] for r in rows])),xmm=np.array([r[4] for r in rows],float))
for k in ('lxi','logM','z','ldd'): cov[k]=cov[k]-np.average(cov[k])   # centre
def design(names,idx=None):
    cols=[np.ones(len(y))]+[cov[n] for n in names]; Xm=np.column_stack(cols); return Xm if idx is None else Xm[idx]
def fit(Xm,yy,ss):
    from scipy.optimize import minimize_scalar
    def nll(lt,ret=False):
        t2=np.exp(2*lt); w=1/(ss**2+t2); A=Xm.T@(Xm*w[:,None]); b=np.linalg.solve(A,Xm.T@(w*yy)); r=yy-Xm@b
        v=0.5*np.sum(np.log(ss**2+t2)+w*r**2)
        return (v,b,np.sqrt(np.diag(np.linalg.inv(A))),np.exp(lt)) if ret else v
    o=minimize_scalar(nll,bounds=(-6,1),method='bounded'); return nll(o.x,True)
models=[("constant",[]),("logM",['logM']),("z",['z']),("ldd",['ldd']),("log xi",['lxi']),("instrument",['xmm']),
        ("logM+z",['logM','z']),("logM+ldd",['logM','ldd']),("log xi+ldd",['lxi','ldd']),("logM+z+ldd",['logM','z','ldd']),("logM+ldd+instr",['logM','ldd','xmm']),("log xi+ldd+instr",['lxi','ldd','xmm'])]
samples=sorted(set(samp))
res=[]
print("\n%-20s k   tau    AIC   dAIC  LOSO_score(sum)  dLOSO  | coefficients"%"model")
for lab,nm in models:
    Xm=design(nm); v,b,se,tau=fit(Xm,y,s); aic=2*v+2*(len(nm)+2)
    loso=0
    for sm in samples:
        te=samp==sm; tr=~te
        vt,bt_,_,tt=fit(Xm[tr],y[tr],s[tr]); mu=Xm[te]@bt_; var=s[te]**2+tt**2
        loso+=np.sum(-0.5*np.log(2*np.pi*var)-0.5*(y[te]-mu)**2/var)
    res.append((lab,nm,aic,loso,b,se,tau))
a0=res[0][2]; l0=res[0][3]
for lab,nm,aic,loso,b,se,tau in res:
    print("%-20s %d  %.3f  %7.1f %+6.1f   %8.1f      %+6.1f | "%(lab,len(nm),tau,aic,aic-a0,loso,loso-l0)+"; ".join("%s=%+.3f±%.3f"%(n,bb,ee) for n,bb,ee in zip(['mu']+nm,b,se)))
# hold-out XXL only (train on all others)
print("\nTrain on non-XXL, predict XXL (162): mean log score per system")
tr=samp!="XXL C1"; te=~tr
for lab,nm in models:
    Xm=design(nm); vt,bt_,_,tt=fit(Xm[tr],y[tr],s[tr]); mu=Xm[te]@bt_; var=s[te]**2+tt**2
    sc=np.mean(-0.5*np.log(2*np.pi*var)-0.5*(y[te]-mu)**2/var); print("  %-20s %.3f   pred mean p=%.2f  (XXL median p_i=%.2f)"%(lab,sc,np.mean(mu),np.median(y[te])))
# ---- XXL only extras (declared: Lx dropped, circular with Mgas through r~0.98 covariance)
raw={}
for ln in open('xxl365gc_vizier.txt'):
    f=ln.rstrip('\n').split('\t')
    if len(f)>=19 and f[0].strip().isdigit():
        try: ng=float(f[5])
        except: ng=np.nan
        raw[int(f[1])]=ng
sel=np.where(samp=="XXL C1")[0]
ng=np.array([raw[int(rows[i][0][5:])] for i in sel]); ok=np.isfinite(ng)&(ng>0)
yy=y[sel][ok]; ss=s[sel][ok]
lng=np.log10(ng[ok]); lM=cov['logM'][sel][ok]; zz=cov['z'][sel][ok]; ld=cov['ldd'][sel][ok]
lng_res=lng-np.polyval(np.polyfit(lM,lng,1),lM)   # richness at fixed mass
print("\nXXL-only (N=%d, single instrument/method): hierarchical fits"%ok.sum())
def show(lab,cols):
    Xm=np.column_stack([np.ones(len(yy))]+cols); v,b,se,tau=fit(Xm,yy,ss); print("  %-26s tau=%.3f AIC=%.1f | "%(lab,tau,2*v+2*(len(cols)+2))+"; ".join("%+.3f±%.3f"%(bb,ee) for bb,ee in zip(b,se)))
show("constant",[]); show("ldd",[ld]); show("logM",[lM]); show("z",[zz]); show("ldd+logM",[ld,lM]); show("ldd+z",[ld,zz]); show("ldd+logM+z",[ld,lM,zz])
show("richness|M",[lng_res]); show("ldd+richness|M",[ld,lng_res])
print("  spread of ldd within XXL: %.2f..%.2f (log10); corr(ldd,logM)=%.2f corr(ldd,z)=%.2f"%(ld.min()+np.average(cov['ldd']),ld.max()+np.average(cov['ldd']),np.corrcoef(ld,lM)[0,1],np.corrcoef(ld,zz)[0,1]))
