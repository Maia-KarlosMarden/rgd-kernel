exec(open('ref41.py').read().split('if __name__=="__main__":')[0])
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
def fit(rows,cols,reml=False):
    y=np.array([r[6] for r in rows]); s=np.array([r[7] for r in rows]); X=np.column_stack([np.ones(len(rows))]+[np.array([r[c] for r in rows],float) for c in cols])
    def comp(lt):
        t2=np.exp(2*lt); w=1/(s**2+t2); A=X.T@(X*w[:,None]); b=np.linalg.solve(A,X.T@(w*y)); r=y-X@b
        return t2,w,A,b,r
    def nll(lt):
        t2,w,A,b,r=comp(lt); v=0.5*np.sum(np.log(s**2+t2)+w*r**2)
        return v+(0.5*np.log(np.linalg.det(A)) if reml else 0)
    o=minimize_scalar(nll,bounds=(-7,1),method="bounded"); t2,w,A,b,r=comp(o.x); return dict(tau=np.sqrt(t2),beta=b,se=np.sqrt(np.diag(np.linalg.inv(A))),nll=o.fun,n=len(rows),chi=np.sum(w*r**2))
def tau_profile(rows,cols):
    best=fit(rows,cols); ts=np.linspace(0,0.4,161); 
    y=np.array([r[6] for r in rows]); s=np.array([r[7] for r in rows]); X=np.column_stack([np.ones(len(rows))]+[np.array([r[c] for r in rows],float) for c in cols])
    out=[]
    for t in ts:
        w=1/(s**2+t**2); A=X.T@(X*w[:,None]); b=np.linalg.solve(A,X.T@(w*y)); r=y-X@b; out.append(0.5*np.sum(np.log(s**2+t**2)+w*r**2))
    out=np.array(out); ok=ts[out<out.min()+0.5*3.84/1]  # 95% (chi2_1/2)
    return ts[np.argmin(out)],ok.min(),ok.max()
def run(rows,tag,verbose=True):
    print("\n=== %s  (N=%d)"%(tag,len(rows)))
    base=fit(rows,[]); t0,lo,hi=tau_profile(rows,[])
    print(" intercept only: mu=%.3f±%.3f  tau=%.3f  (profile 95%% [%.3f, %.3f])  chi2/N of the mean model (tau=0) -> see below"%(base['beta'][0],base['se'][0],base['tau'],lo,hi))
    y=np.array([r[6] for r in rows]); s=np.array([r[7] for r in rows]); w=1/s**2; mu=np.sum(w*y)/np.sum(w); print("   homogeneity (tau=0): chi2=%.0f for %d dof (p=%.1e)"%(np.sum(w*(y-mu)**2),len(y)-1,chi2.sf(np.sum(w*(y-mu)**2),len(y)-1)))
    for c in (3,2,4,5):
        f=fit(rows,[c]); lr=2*(base['nll']-f['nll']); print("   +%-20s beta=%+.3f±%.3f  tau %.3f -> %.3f   LRT=%.1f (p=%.3g)"%(names[c],f['beta'][1],f['se'][1],base['tau'],f['tau'],lr,chi2.sf(max(lr,0),1)))
    full=fit(rows,[3,2,4,5]); lr=2*(base['nll']-full['nll']); print("   all four: "+"; ".join("%s=%+.3f±%.3f"%(l,b,e) for l,b,e in zip(["mu","logM","z","XMM","logDelta"],full['beta'],full['se']))+"  tau=%.3f  LRT=%.1f (p=%.3g, 4 dof)"%(full['tau'],lr,chi2.sf(max(lr,0),4)))
    return base,full
rows=build(); run(rows,"A. all measured-f_gas samples, alpha=1, sigma_f=10% (Allen: own errors), sigma_M=10%")
# bootstrap over unique clusters for the full model
def cluster_boot(rows,B=1000,seed=2):
    rng=np.random.default_rng(seed); ids=sorted(set(norm(r[0]) for r in rows)); groups={i:[r for r in rows if norm(r[0])==i] for i in ids}
    out=[]
    for _ in range(B):
        pick=rng.choice(len(ids),len(ids)); rr=[x for k in pick for x in groups[ids[k]]]
        try: f=fit(rr,[3,2,4,5]); b=fit(rr,[]); out.append(np.r_[f['beta'],f['tau'],b['tau'],b['beta'][0]])
        except Exception: pass
    return np.array(out)
bs=cluster_boot(rows); print("\n cluster bootstrap (N unique clusters = %d, 1000 resamples), 95%% intervals:"%len(set(norm(r[0]) for r in rows)))
for j,l in enumerate(["mu","beta_logM (per dex)","beta_z","beta_XMM","beta_logDelta","tau (full)","tau (intercept-only)","mu (intercept-only)"]): print("   %-22s %+.3f  [%+.3f, %+.3f]"%(l,np.median(bs[:,j]),*np.percentile(bs[:,j],[2.5,97.5])))
# variants
run(build(sfmult=0.5),"B. smaller measurement errors (sigma_f x0.5)")
run(build(sfmult=1.5),"C. larger measurement errors (sigma_f x1.5)")
run([r for r in rows if not r[8]],"D. drop the 3 X-COP mergers")
seen={}; ded=[]
for pref in ("Ettori r200","Vikhlinin","Allen"):
    pass
order={"Ettori r200":0,"Vikhlinin":1,"Allen":2,"Groups-L(XMM)":3,"Groups-S(Chandra)":3}
best={}
for r in rows:
    k=norm(r[0]); 
    if k not in best or order[r[1]]<order[best[k][1]]: best[k]=r
run(list(best.values()),"E. one measurement per cluster (priority Ettori > Vikhlinin > Allen > groups)")
for a in (0.5,1.5): run(build(alpha=a),"F. g_dagger ~ E^%.1f (degeneracy of the z trend with alpha)"%a)
run([r for r in rows if r[1] in("Vikhlinin","Allen","Ettori r200")],"G. clusters only (no groups)")

def run2(rows,tag,cs):
    print("\n=== %s  (N=%d)"%(tag,len(rows))); base=fit(rows,[])
    print(" intercept only: mu=%.3f±%.3f tau=%.3f"%(base['beta'][0],base['se'][0],base['tau']))
    for c in cs:
        f=fit(rows,[c]); lr=2*(base['nll']-f['nll']); print("   +%-20s beta=%+.3f±%.3f tau %.3f -> %.3f LRT=%.1f (p=%.3g)"%(names[c],f['beta'][1],f['se'][1],base['tau'],f['tau'],lr,chi2.sf(max(lr,0),1)))
run2([r for r in rows if r[1].startswith("Groups")],"H. groups only (Lovisari XMM vs Sun Chandra)",(3,2,4))
# ---- I. hierarchical model with the switch: p_i = p0 / (1+(xi_c/xi_i)^m) + eps
def switch_fit(rows,extra=[]):
    y=np.array([r[6] for r in rows]); s=np.array([r[7] for r in rows]); xi_=np.array([r[9] for r in rows])
    Xe=np.column_stack([np.array([r[c] for r in rows],float) for c in extra]) if extra else np.zeros((len(rows),0))
    best=None
    for m in (1,2,3,4,6,8,16,32):
        for xc in np.exp(np.linspace(np.log(.04),np.log(.5),60)):
            th=1/(1+(xc/xi_)**m); X=np.column_stack([th,Xe])
            def nll(lt):
                t2=np.exp(2*lt); w=1/(s**2+t2); A=X.T@(X*w[:,None]); b=np.linalg.solve(A,X.T@(w*y)); r=y-X@b; return 0.5*np.sum(np.log(s**2+t2)+w*r**2),b,t2
            o=minimize_scalar(lambda lt: nll(lt)[0],bounds=(-7,1),method="bounded"); v,b,t2=nll(o.x)
            if best is None or v<best[0]: best=(v,m,xc,b,np.sqrt(t2))
    return best
print("\n=== I. hierarchical model with the xi switch (all %d data): p_i = p0*Theta(xi_i;xi_c,m)+eps, eps~N(0,tau^2+sigma_i^2)"%len(rows))
base=fit(rows,[]); 
sw=switch_fit(rows); print(" constant p      : -lnL %.1f  tau=%.3f  (2 par)"%(base['nll'],base['tau']))
print(" switch (m, xi_c, p0): -lnL %.1f  m=%d xi_c=%.3f p0=%.3f tau=%.3f  (4 par)  LRT vs constant = %.1f"%(sw[0],sw[1],sw[2],sw[3][0],sw[4],2*(base['nll']-sw[0])))
for extra in ([3],[2],[4],[5],[3,2,4,5]):
    f=switch_fit(rows,extra); print("   switch + %-32s: tau=%.3f  betas %s  (m=%d xi_c=%.3f p0=%.3f) LRT vs switch alone=%.1f"%(",".join(names[c] for c in extra),f[4],np.round(f[3][1:],3),f[1],f[2],f[3][0],2*(sw[0]-f[0])))
cl=[r for r in rows if r[1] in("Vikhlinin","Allen","Ettori r200")]
swc=switch_fit(cl); bc=fit(cl,[]); print(" clusters only: constant tau=%.3f (-lnL %.1f); switch tau=%.3f (-lnL %.1f, m=%d xi_c=%.2f)"%(bc['tau'],bc['nll'],swc[4],swc[0],swc[1],swc[2]))

