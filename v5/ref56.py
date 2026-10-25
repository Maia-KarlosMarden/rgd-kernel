"""PRE-REGISTERED (written before running).  Paired comparison of RGD v6 (x_eff=x(xi/y)^theta, q=1) against RGD v5 (s=(Delta/Delta500)^(0.5*Theta)), both with the SAME Theta(xi_c=.115,m=4),
no parameter refit on these data.  Systems (non-circular; X-COP constant-f excluded): Vikhlinin r2500 (10), Ettori r200 (12), Allen r2500 (37+), groups r2500 (43), XXL-365 C1 at 500 kpc (162).
Per system: r_i = ln f_obs - ln f_pred, variance sf_i^2 + tau^2 (tau refit per model on the pooled set).  Statistic: DeltaLS = LS(v6) - LS(v5) (sum of log scores; >0 favours v6).
Uncertainty: bootstrap of systems stratified by sample (4000).  Also: win fraction (|r6|<|r5|), per-sample table, equal-sample-weight mean, 1-parameter versions (best p vs best q),
and variants: Theta=1; XXL NFW concentration x0.7 / x1.3.  v6 was chosen AFTER seeing ref52 (ratio family), so this is a descriptive comparison, not a blind test (stated up front)."""
import io,contextlib,sys,numpy as np,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import ref52 as T
R=T.R; X=T.Q.X; Ms=R.Ms; Mpc=R.Mpc
from scipy.optimize import minimize_scalar
XIC,MW=0.115,4
def theta(xi): return 1/(1+(XIC/xi)**MW)
def prep(recs):
    out=[]
    for d in recs:
        gd=T.gd_of(d); xi=T.Gc*d['M5']/d['R5']**2/gd; dd=(d['M']/d['M5'])/(d['r']/d['R5'])**3; gr=(d['M']/d['M5'])/(d['r']/d['R5'])**2
        out.append(dict(samp=d['samp'],xi=xi,dd=dd,gr=gr,rec=d))
    return out
def resid(P,law,p=0.5,q=1.0,thetaone=False):
    r=np.empty(len(P))
    for i,e in enumerate(P):
        th=1.0 if thetaone else theta(e['xi'])
        s=e['dd']**(p*th) if law=='v5' else e['gr']**(q*th)
        fp=T.fpred_s(e['rec'],s); r[i]=np.log(e['rec']['fb'])-np.log(fp) if np.isfinite(fp) and fp>0 else np.nan
    return r
def tau_fit(r,sf):
    o=minimize_scalar(lambda lt: 0.5*np.sum(np.log(2*np.pi*(sf**2+np.exp(2*lt)))+r**2/(sf**2+np.exp(2*lt))),bounds=(-6,1),method='bounded'); return np.exp(o.x),o.fun
def ls_i(r,sf,tau): v=sf**2+tau**2; return -0.5*(np.log(2*np.pi*v)+r**2/v)
def compare(P,label,thetaone=False,boot=4000):
    sf=np.array([e['rec']['sf'] for e in P]); samp=np.array([e['samp'] for e in P])
    r5=resid(P,'v5',thetaone=thetaone); r6=resid(P,'v6',thetaone=thetaone); ok=np.isfinite(r5)&np.isfinite(r6)
    r5,r6,sf,samp=r5[ok],r6[ok],sf[ok],samp[ok]
    t5,_=tau_fit(r5,sf); t6,_=tau_fit(r6,sf); l5=ls_i(r5,sf,t5); l6=ls_i(r6,sf,t6); d=l6-l5
    rng=np.random.default_rng(7); idxs={s:np.where(samp==s)[0] for s in set(samp)}; bs=[]
    for _ in range(boot):
        ii=np.concatenate([rng.choice(ix,len(ix)) for ix in idxs.values()])
        a5,_=tau_fit(r5[ii],sf[ii]); a6,_=tau_fit(r6[ii],sf[ii]); bs.append(ls_i(r6[ii],sf[ii],a6).sum()-ls_i(r5[ii],sf[ii],a5).sum())
    bs=np.array(bs); win=np.mean(np.abs(r6)<np.abs(r5))
    wm=np.mean([d[samp==s].mean() for s in idxs])
    print("\n== %s  (N=%d)  tau5=%.3f tau6=%.3f | rms ln: v5 %.3f  v6 %.3f"%(label,len(r5),t5,t6,np.sqrt(np.mean(r5**2)),np.sqrt(np.mean(r6**2))))
    print("   DeltaLS(v6-v5)=%+.1f  [95%% boot %+.1f, %+.1f]   P(DeltaLS>0)=%.3f   win fraction v6=%.2f   equal-sample-weight mean per-system dLS=%+.3f"%(d.sum(),*np.percentile(bs,[2.5,97.5]),np.mean(bs>0),win,wm))
    print("   %-18s N   rms5   rms6   mean r5  mean r6   dLS/system"%"sample")
    for s in sorted(idxs): 
        m=samp==s; print("   %-18s %3d  %.3f  %.3f  %+.3f   %+.3f    %+.3f"%(s,m.sum(),np.sqrt(np.mean(r5[m]**2)),np.sqrt(np.mean(r6[m]**2)),r5[m].mean(),r6[m].mean(),d[m].mean()))
    return d.sum(),np.percentile(bs,[2.5,97.5])
P0=prep(T.D)
P0=[e for e in P0 if e['samp']!="__"]
compare(P0,"ALL, archived Theta")
compare([e for e in P0 if e['samp']!="XXL C1"],"WITHOUT XXL (Vikhlinin, Ettori, Allen, groups)")
compare([e for e in P0 if e['samp']=="XXL C1"],"XXL C1 only")
compare(P0,"ALL, Theta=1",thetaone=True)
# concentration variants for XXL
base=X.c200_DM14
def rebuild(fac):
    X.c200_DM14=lambda M,z,h=.7,fac=fac: fac*base(M,z,h)
    new=[]
    for d in X.rows:
        if not (np.isfinite(d['Mg']) and np.isfinite(d['r500']) and d['Mg']>0 and d['cls']==1): continue
        M5=X.M500_from_r(d['r500'],d['z']); r5=d['r500']*1e-3
        try: Mr=X.M_in(M5,r5,d['z'])
        except Exception: continue
        st=R.Mstar_tot(M5)*Ms; Mrk=Mr*Ms
        new.append(dict(name="XLSSC%d"%d['id'],samp="XXL C1",z=d['z'],instr=1,M=Mrk,r=0.5*Mpc,M5=M5*Ms,R5=r5*Mpc,dd=(Mr/M5)/(0.5/r5)**3,
                        fb=d['Mg']*Ms/Mrk+st/Mrk,sf=np.hypot(d['eMg']*1e11/d['Mg'],0.15),H0=R.H0))
    X.c200_DM14=base; return new
for fac in (0.7,1.3):
    P=[e for e in P0 if e['samp']!="XXL C1"]+prep(rebuild(fac)); compare(P,"ALL, XXL concentration x%.1f"%fac,boot=1500)
# one-parameter versions: best p (v5) vs best q (v6) on ALL
sfA=np.array([e['rec']['sf'] for e in P0])
def nllp(par,law):
    r=resid(P0,law,p=par,q=par); ok=np.isfinite(r)
    return tau_fit(r[ok],sfA[ok])[1]
from scipy.optimize import minimize_scalar as ms
b5=ms(lambda p: nllp(p,'v5'),bounds=(0.2,1.5),method='bounded'); b6=ms(lambda q: nllp(q,'v6'),bounds=(0.2,2.0),method='bounded')
print("\n1-parameter versions on ALL: v5 best p=%.3f nll=%.1f  AIC=%.1f | v6 best q=%.3f nll=%.1f AIC=%.1f | dAIC(v6-v5)=%+.1f"%(b5.x,b5.fun,2*b5.fun+4,b6.x,b6.fun,2*b6.fun+4,2*(b6.fun-b5.fun)))
# --- declared addition: which samples drive the 1-parameter advantage?  (fits on subsets)
def one_param(P,label):
    sfx=np.array([e['rec']['sf'] for e in P])
    def nl(par,law):
        r=resid(P,law,p=par,q=par); ok=np.isfinite(r); return tau_fit(r[ok],sfx[ok])[1]
    a=ms(lambda p: nl(p,'v5'),bounds=(0.2,1.5),method='bounded'); b=ms(lambda q: nl(q,'v6'),bounds=(0.2,2.0),method='bounded')
    print("1-param %-26s v5 p=%.3f AIC=%.1f | v6 q=%.3f AIC=%.1f | dAIC(v6-v5)=%+.1f"%(label,a.x,2*a.fun+4,b.x,2*b.fun+4,2*(b.fun-a.fun)))
one_param([e for e in P0 if e['samp']!="XXL C1"],"WITHOUT XXL")
one_param([e for e in P0 if e['samp']=="XXL C1"],"XXL only")
