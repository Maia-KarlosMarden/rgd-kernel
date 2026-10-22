exec(open('ref41b.py').read().split("rows=build(); run(rows")[0])
from scipy.stats import spearmanr,kendalltau
def zfit(rr):
    f=fit(rr,[2]); b=fit(rr,[]); return f['beta'][1],f['se'][1],f['tau'],2*(b['nll']-f['nll'])
CL=("Vikhlinin","Allen","Ettori r200")
print("RAW CORRELATION of p_i with z (alpha=1)")
rows=build()
for lab,sel in (("clusters (59)",lambda r:r[1] in CL),("Allen only",lambda r:r[1]=="Allen"),("Allen z>0.3",lambda r:r[1]=="Allen" and r[2]>0.3),("all 102",lambda r:True)):
    rr=[r for r in rows if sel(r)]; z_=np.array([r[2] for r in rr]); p_=np.array([r[6] for r in rr])
    sp=spearmanr(z_,p_); kt=kendalltau(z_,p_)
    rng=np.random.default_rng(1); bs=[spearmanr(z_[i],p_[i])[0] for i in (rng.integers(0,len(rr),len(rr)) for _ in range(2000))]
    print(" %-14s N=%3d z %.2f-%.2f  Spearman rho=%+.2f (p=%.3g) [boot 95%% %+.2f,%+.2f] ; Kendall tau=%+.2f (p=%.3g)"%(lab,len(rr),z_.min(),z_.max(),sp[0],sp[1],*np.percentile(bs,[2.5,97.5]),kt[0],kt[1]))
print("\nHIERARCHICAL slope of p with z (per unit z), clusters only, as a function of alpha in g_dagger ~ E(z)^alpha")
print("  alpha | N  | beta_z ± se   | tau   | LRT(beta_z=0) | mu")
res=[]
for a in (0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.5):
    rr=[r for r in build(alpha=a) if r[1] in CL]; b,se,t,lr=zfit(rr); mu=fit(rr,[2])['beta'][0]; res.append((a,b,se)); print("  %.2f  | %2d | %+.3f±%.3f | %.3f | %5.1f | %.3f"%(a,len(rr),b,se,t,lr,mu))
res=np.array(res); 
for i in range(len(res)-1):
    if res[i,1]*res[i+1,1]<0:
        a0=res[i,0]+(res[i+1,0]-res[i,0])*(-res[i,1])/(res[i+1,1]-res[i,1]); print("  beta_z=0 at alpha = %.2f"%a0)
print("\nSame, Allen only (single instrument, z=0.06-1.1, r2500):")
for a in (0.0,0.5,1.0,1.5,2.0):
    rr=[r for r in build(alpha=a) if r[1]=="Allen"]; b,se,t,lr=zfit(rr); print("  alpha=%.1f N=%d beta_z=%+.3f±%.3f tau=%.3f LRT=%.1f"%(a,len(rr),b,se,t,lr))
print("\nDoes the slope persist with covariates (clusters, alpha=1)? logM + z:")
rr=[r for r in rows if r[1] in CL]; f=fit(rr,[3,2]); print("  logM %+.3f±%.3f ; z %+.3f±%.3f ; tau %.3f"%(f['beta'][1],f['se'][1],f['beta'][2],f['se'][2],f['tau']))
# bootstrap over clusters for the z slope at alpha=1, and a high-z/low-z split
rng=np.random.default_rng(5); ids=sorted(set(norm(r[0]) for r in rr)); g={i:[r for r in rr if norm(r[0])==i] for i in ids}
bs=[]
for _ in range(1000):
    pick=rng.choice(len(ids),len(ids)); q=[x for k in pick for x in g[ids[k]]]
    try: bs.append(fit(q,[2])['beta'][1])
    except Exception: pass
print("  cluster bootstrap beta_z (alpha=1): median %+.3f [%+.3f, %+.3f]"%(np.median(bs),*np.percentile(bs,[2.5,97.5])))
lo=[r for r in rr if r[2]<0.15]; hi=[r for r in rr if r[2]>=0.15]
for lab,q in (("z<0.15",lo),("z>=0.15",hi)):
    b=fit(q,[]); print("  %-8s N=%2d mu=%.3f±%.3f tau=%.3f  median z %.2f"%(lab,len(q),b['beta'][0],b['se'][0],b['tau'],np.median([r[2] for r in q])))
