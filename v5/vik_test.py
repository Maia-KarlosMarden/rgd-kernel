import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
from scipy.optimize import brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('stars_test.py').read().split('print("\\nALLEN')[0])
# Vikhlinin+2006 Table 3 (+Table 1 z); h=0.71, Om=0.3 ; M in 1e14 Msun, r500 kpc, fgas cumulative (h^-3/2)
V=[("A133",.0569,998,1.24,3.14,.065,.084),("A262",.0162,663,.34,None,.067,None),("A383",.1883,956,1.68,3.10,.090,.122),("A478",.0881,1359,4.23,7.83,.096,.118),
("A907",.1603,1117,2.30,4.71,.090,.121),("A1413",.1429,1339,3.08,7.78,.092,.104),("A1795",.0622,1283,2.66,6.57,.088,.098),("A1991",.0592,753,.60,1.28,.068,.099),
("A2029",.0779,1380,4.41,8.29,.090,.121),("A2390",.2302,1448,3.50,10.88,.127,.138),("RXJ1159",.0810,695,.31,None,.042,None),("MKW4",.0199,635,.30,.74,.044,.063)]
Gc=6.674e-11; cc=3e8; Ms=1.989e30; kpc=3.0857e19; Mpc=3.0857e22; lam_=.999; A_=.4068
H=lambda z: 71e3/Mpc*np.sqrt(.3*(1+z)**3+.7)
def solve(m,r,xi,gd,s,n=None):
    g=Gc*m/r**2; K=lambda f: 1+lam_*np.exp(-f*g/gd/s)/(-np.expm1(-f*g/gd/s))**(.5+A_*xi); return brentq(lambda f: f*K(f)-1,1e-4,5.0)
rows_=[]
for n,z,r5,M25,M5,f25,f5 in V:
    rc_=3*H(z)**2/(8*np.pi*Gc); r500=r5*kpc; M500d=500*rc_*4*np.pi/3*r500**3
    M25_=M25*1e14*Ms; r25=(M25_/(2500*rc_*4*np.pi/3))**(1/3)
    M500=(M5*1e14*Ms) if M5 else M500d
    xi=Gc*M500/r500**2/(cc*H(z)/(2*np.pi)); rows_.append(dict(n=n,z=z,r500=r500,M500=M500,M25=M25_,r25=r25,f25=f25,f5=f5,xi=xi,gd=cc*H(z)/(2*np.pi),chk=(M5*1e14*Ms/M500d) if M5 else np.nan,w=r25/r500))
print("consistency M500/[500 rho_c 4pi/3 r500^3] (reported M500): "+", ".join("%s %.2f"%(r['n'],r['chk']) for r in rows_ if np.isfinite(r['chk'])))
print("clusters with both radii:",sum(1 for r in rows_ if r['f5']),"; xi range %.2f-%.2f ; r2500/r500 range %.2f-%.2f ; M500 %.1e-%.1e Msun"%(min(r['xi'] for r in rows_),max(r['xi'] for r in rows_),min(r['w'] for r in rows_),max(r['w'] for r in rows_),min(r['M500'] for r in rows_)/Ms,max(r['M500'] for r in rows_)/Ms))
PHI=.7
models={"p=0.50":lambda x:.5,"p=0.46":lambda x:.46,"p(xi)=0.46(xi/0.4)^0.37":lambda x:.46*(x/.4)**.37}
ST500=lambda r: Mstar_tot(r['M500']/Ms)*Ms
print("\n%-8s %5s %5s | M*/Mgas(500) (2500)| f_bar obs: r500  r2500 ratio | pred r500 | pred r2500 [p=.5 / p=.46 / p(xi)] | obs/pred at 2500 [p=.5 / p=.46 / p(xi)]"%("cluster","xi","z"))
res={k:[] for k in models}; res500=[]; lm=[]; ratio_obs=[]; ratio_pred={k:[] for k in models}
for r in rows_:
    st25=PHI*ST500(r)/(r['f25']*r['M25']); fb25=r['f25']*(1+st25)
    out=[];pr25=[]
    for k,pf_ in models.items():
        p=pf_(r['xi']); f25p=solve(r['M25'],r['r25'],r['xi'],r['gd'],5.0**p); pr25.append(f25p); res[k].append(np.log(fb25/f25p))
    lm.append(np.log10(r['M25']/Ms))
    if r['f5']:
        st5=ST500(r)/(r['f5']*r['M500']); fb5=r['f5']*(1+st5); f5p=solve(r['M500'],r['r500'],r['xi'],r['gd'],1.0); res500.append(np.log(fb5/f5p)); ratio_obs.append(fb25/fb5)
        for i,k in enumerate(models): ratio_pred[k].append(pr25[i]/f5p)
        print("%-8s %5.2f %5.3f | %.2f (%.2f) | %.3f %.3f %.2f | %.3f | %.3f %.3f %.3f | %.2f %.2f %.2f"%(r['n'],r['xi'],r['z'],st5,st25,fb5,fb25,fb25/fb5,f5p,*pr25,*[fb25/x for x in pr25]))
    else:
        print("%-8s %5.2f %5.3f | (no r500 f_gas) %.2f | f_bar2500 obs %.3f | pred %.3f %.3f %.3f"%(r['n'],r['xi'],r['z'],st25,fb25,*pr25))
lm=np.array(lm); print("\nAt r500 (w=1; p irrelevant): obs/pred median %.2f, rms(ln) %.3f (N=%d)"%(np.exp(np.median(res500)),np.sqrt(np.mean(np.square(res500))),len(res500)))
print("At r2500:")
for k in models:
    a=np.array(res[k]); print("  %-26s obs/pred median %.2f, rms(ln) %.3f ; trend of ln(obs/pred) with log M2500: %+.2f"%(k,np.exp(np.median(a)),np.sqrt(np.mean(a**2)),np.polyfit(lm,a,1)[0]))
ro=np.array(ratio_obs); print("Within-cluster ratio f_bar(r2500)/f_bar(r500), 10 clusters: observed median %.2f ; predicted: "%np.median(ro)+" ; ".join("%s %.2f"%(k,np.median(v)) for k,v in ratio_pred.items()))
lm10=np.array([np.log10(r['M25']/Ms) for r in rows_ if r['f5']])
print("  trend of the ratio with log M2500: observed %+.2f ; predicted: "%np.polyfit(lm10,np.log(ro),1)[0]+" ; ".join("%s %+.2f"%(k,np.polyfit(lm10,np.log(v),1)[0]) for k,v in ratio_pred.items()))
