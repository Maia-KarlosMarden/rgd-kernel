src=open('ref41b.py').read().split("rows=build(); run(rows")[0]
exec(src)
exec("def switch_fit"+open('ref41b.py').read().split("def switch_fit")[1].split("print(\"\\n=== I.")[0])
rows=build()
sw=switch_fit(rows); m0,xc0=sw[1],sw[2]
def tau_with_switch(rr,m,xc,extra=[]):
    y=np.array([r[6] for r in rr]); s=np.array([r[7] for r in rr]); th=1/(1+(xc/np.array([r[9] for r in rr]))**m)
    Xe=np.column_stack([np.array([r[c] for r in rr],float) for c in extra]) if extra else np.zeros((len(rr),0)); X=np.column_stack([th,Xe])
    def nll(lt):
        t2=np.exp(2*lt); w=1/(s**2+t2); A=X.T@(X*w[:,None]); b=np.linalg.solve(A,X.T@(w*y)); r=y-X@b; return 0.5*np.sum(np.log(s**2+t2)+w*r**2),b,t2,np.linalg.inv(A)
    o=minimize_scalar(lambda lt: nll(lt)[0],bounds=(-7,1),method="bounded"); v,b,t2,cov=nll(o.x); return np.sqrt(t2),b,np.sqrt(np.diag(cov))
for lab,sel in (("all",lambda r:True),("clusters",lambda r:r[1] in("Vikhlinin","Allen","Ettori r200")),("groups",lambda r:r[1].startswith("Groups")),("groups xi>=0.2",lambda r:r[1].startswith("Groups") and r[9]>=0.2)):
    rr=[r for r in rows if sel(r)]; t,b,se=tau_with_switch(rr,m0,xc0); print("%-16s N=%3d  with switch (m=%d, xi_c=%.3f): p0=%.3f±%.3f  tau=%.3f"%(lab,len(rr),m0,xc0,b[0],se[0],t))
print("\nresidual covariate effects after the switch, one at a time (SE, GLS):")
for c in (3,2,4,5):
    t,b,se=tau_with_switch(rows,m0,xc0,[c]); print("   %-20s beta=%+.3f±%.3f  tau=%.3f"%(names[c],b[1],se[1],t))
rng=np.random.default_rng(7); ids=sorted(set(norm(r[0]) for r in rows)); g={i:[r for r in rows if norm(r[0])==i] for i in ids}
out=[]
for _ in range(300):
    pick=rng.choice(len(ids),len(ids)); rr=[x for k in pick for x in g[ids[k]]]
    f=switch_fit(rr); base=fit(rr,[]); out.append([f[3][0],f[2],f[1],f[4],base['tau'],2*(base['nll']-f[0])]); 
out=np.array(out)
print("\ncluster bootstrap (300): p0 %.3f [%.3f,%.3f] ; xi_c median %.3f [%.3f,%.3f] ; tau(switch) %.3f [%.3f,%.3f] ; tau(const) %.3f [%.3f,%.3f] ; LRT switch vs const: median %.1f, 2.5%% quantile %.1f"%(np.median(out[:,0]),*np.percentile(out[:,0],[2.5,97.5]),np.median(out[:,1]),*np.percentile(out[:,1],[2.5,97.5]),np.median(out[:,3]),*np.percentile(out[:,3],[2.5,97.5]),np.median(out[:,4]),*np.percentile(out[:,4],[2.5,97.5]),np.median(out[:,5]),np.percentile(out[:,5],2.5)))
