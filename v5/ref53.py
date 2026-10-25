"""PRE-REGISTERED (written before running). Stress test of the parameter-free law s = g(r)/g500 on X-COP (13 clusters x 3 radii: 0.5,1,1.5 Mpc, resolved profiles).
With s=g/g500 the kernel argument is f*xi at every radius -> f_bar(r) would be radius-independent inside each cluster. X-COP data show whether that is acceptable.
Compare rms of log10(f_obs/f_pred) for: none, v4, s=dd^0.5 (archived), s=g/g500 (q=1), s=(g/g500)^q (q grid).  Same data as ref31 section A."""
import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
buf=io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(open('stars_test.py').read())
    cd={}; exec(open('kernel_scenarios.py').read().split("# ---------- check NFW R500")[0],cd)
D0=cd['build'](radii_mpc=[.5,1,1.5]); lam=cd.get('lam',.999); A=cd.get('A',.4068)
xi=D0['xi']; w=D0['w']; fb=D0['fb']; g=D0['g']; gd=D0['gd']; Mr=D0['Mr']; M5=D0['M5']
dd=(Mr/M5)/w**3; gr=(Mr/M5)/w**2
x=fb*g/gd
K=lambda n,s: 1+lam*np.exp(-x/s)/(-np.expm1(-x/s))**n
rms=lambda a: float(np.sqrt(np.nanmean(np.square(a))))
res=lambda n,s: np.log10(K(n,s)*fb)
n5=.5+A*xi
print("N points",len(xi))
print("none            rms %.3f dex"%rms(res(n5,1.0)))
print("dd^0.5 (arch.)  rms %.3f dex"%rms(res(n5,dd**.5)))
for p in (.45,.6,.8): print("dd^%.2f          rms %.3f dex"%(p,rms(res(n5,dd**p))))
for q in (0.25,0.5,0.75,1.0,1.125):
    print("(g/g500)^%.3f   rms %.3f dex"%(q,rms(res(n5,gr**q))))
# radial behaviour: mean log residual by radius
for lab,s in (("dd^0.5",dd**.5),("g/g500",gr)):
    r=res(n5,s).reshape(-1,3) if len(xi)%3==0 else None
    if r is not None: print("%-8s mean resid at 0.5/1/1.5 Mpc (point order as built): %s"%(lab,np.round(np.nanmean(r,axis=0),3)))
print("observed fb mean by column:",np.round(np.nanmean(fb.reshape(-1,3),axis=0),3),"  gr by column:",np.round(np.nanmean(gr.reshape(-1,3),axis=0),2),"dd:",np.round(np.nanmean(dd.reshape(-1,3),axis=0),2))
