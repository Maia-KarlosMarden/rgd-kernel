"""PRE-REGISTERED (written before running).  Offset-free radial test.  For clusters/groups measured at BOTH r500 and a second radius in the SAME sample
(Vikhlinin: r500,r2500 (10); Ettori: r500,r200 (9); groups Lovisari+Sun: r500,r2500 (43)), the common calibration offset (instrument, normalisation of gas mass) cancels in
   rho_i = ln[ f_obs(r)/f_obs(r500) ] - ln[ f_pred(r)/f_pred(r500) ].
Report mean rho +- s.e.m. and rms per sample for: v5 (p=0.5, 0.55), v6 (q=1, 1.08), and no radial term (s=1), all with the archived Theta, stars Kravtsov (phi=0.7 at r2500 for V/groups, 1 elsewhere, as in the archived prescription).
Ideal model: mean rho = 0.  This test is insensitive to a constant multiplicative offset between samples/instruments."""
import io,contextlib,sys,numpy as np,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import ref52 as T
    import ref55 as H
R=T.R; D=R.D; Ms=R.Ms; Mpc=R.Mpc; Mstar_tot=R.Mstar_tot
E=lambda z: np.sqrt(.3*(1+z)**3+.7)
XIC,MW=0.115,4
th=lambda xi: 1/(1+(XIC/xi)**MW)
# 500-level observed baryon fractions (phi=1)
f500={}
for n,z,r5,M25,M5,f25,f5 in H.V: f500[("Vikhlinin",n)]=f5+Mstar_tot(M5*1e14)/(M5*1e14)
for n,z,M5,M2,f5,f2 in H.Xe: f500[("Ettori r200",n)]=f5+Mstar_tot(M5*1e14)/(M5*1e14)
for r in H.GR_:
    M5=r['M5']/R.Ms; f500[("Groups-L(XMM)" if r['src']=="L" else "Groups-S(Chandra)",r['name'])]=r['f5']+Mstar_tot(M5)/M5
def pred(d,law,par,r500=False):
    gd=T.gd_of(d); xi=T.Gc*d['M5']/d['R5']**2/gd
    if r500: return T.fpred_s(dict(d,M=d['M5'],r=d['R5']),1.0)
    dd=(d['M']/d['M5'])/(d['r']/d['R5'])**3; gr=(d['M']/d['M5'])/(d['r']/d['R5'])**2; t=th(xi)
    s=1.0 if law=='none' else (dd**(par*t) if law=='v5' else gr**(par*t)); return T.fpred_s(d,s)
models=[('none',0),('v5',0.5),('v5',0.55),('v6',1.0),('v6',1.08)]
out={}
for d in D:
    k=(d['samp'],d['name'])
    if k not in f500: continue
    for law,par in models:
        a=pred(d,law,par); b=pred(d,law,par,True)
        if np.isfinite(a) and np.isfinite(b):
            rho=np.log(d['fb']/f500[k])-np.log(a/b); out.setdefault((d['samp'],law,par),[]).append(rho)
print("offset-free ratio test: rho = ln[f_obs(r)/f_obs(r500)] - ln[f_pred(r)/f_pred(r500)]  (0 = perfect)\n")
print("%-18s %-5s %5s   N   mean±sem      rms"%("sample","model","par"))
samples=sorted(set(k[0] for k in out))
for s in samples:
    for law,par in models:
        v=np.array(out[(s,law,par)]); print("%-18s %-5s %5.2f  %3d  %+.3f±%.3f  %.3f"%(s,law,par,len(v),v.mean(),v.std(ddof=1)/np.sqrt(len(v)),np.sqrt(np.mean(v**2))))
    print()
print("pooled (equal weight per cluster):")
for law,par in models:
    v=np.concatenate([out[(s,law,par)] for s in samples]); print("  %-5s %5.2f  N=%d  mean %+.3f±%.3f  rms %.3f"%(law,par,len(v),v.mean(),v.std(ddof=1)/np.sqrt(len(v)),np.sqrt(np.mean(v**2))))
# direct cross-instrument check on shared clusters (r500 baryon fraction)
print("\nShared clusters, same quantity f_b(<r500) from different instruments/methods:")
for n in ("A1795","A2029"):
    a=f500[("Vikhlinin",n)]; b=f500[("Ettori r200",n)]; print("  %-6s Vikhlinin(Chandra) %.3f  Ettori(XMM) %.3f  ratio XMM/Chandra %.2f"%(n,a,b,b/a))
