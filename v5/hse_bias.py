import io,contextlib,numpy as np
from scipy.optimize import minimize, brentq
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('stars_test.py').read())          # Allen arrays: z,r25,fo,M25,M5,xi,lM ; Mstar_tot ; fsolve ; constants
ef=np.array([a[4] for a in ALLEN]); Mgas=fo*M25
def model(th,phi,p,ret=False):
    b0,b1=th; b=np.clip(b0+b1*(lM-14.4),-0.5,0.8); Mt=M25/(1-b)
    Mst=phi*Mstar_tot(M5/Msun/(1-b)**1.0)*Msun        # M*tot from true M500 (approx: M500 scales like M2500)
    fb_true=(Mgas+Mst)/Mt; xit=xi*(1/(1-b))**(1/3)
    pr=np.array([fsolve(m,r,x,c_*Hz(zz)/(2*np.pi),p) for m,r,x,zz in zip(Mt,r25,xit,z)])
    err=ef*(1+Mst/Mgas)/(1-b)
    chi=((fb_true-pr)/(err*(1-b)))**2
    if ret: return fb_true,pr,b,np.sum(chi)
    return np.sum(chi)
print("Allen+2008, 42 clusters. Observed M_HSE = (1-b) M_true ; b(M)=b0+b1*(log M2500 - 14.4). RGD applied to true mass; Mgas measured independently of HSE; M* from Kravtsov relation (phi = fraction inside r2500).")
print("M2500 range: log %.2f-%.2f"%(lM.min(),lM.max()))
for phi in (.5,.7,1.0):
    for p in (.46,.5):
        c0=model([0,0],phi,p); r1=minimize(lambda t: model([t[0],0],phi,p),[0.1],method="Nelder-Mead"); r2=minimize(lambda t: model(t,phi,p),[0.1,0.1],method="Nelder-Mead",options=dict(xatol=1e-3,fatol=1e-3))
        fb,pr,b,c2=model(r2.x,phi,p,ret=True); lo=lM<np.median(lM)
        # equivalent: bias at M2500=2e14, 5e14, 1e15 and the mass-dependent change over the sample
        bf=lambda L: r2.x[0]+r2.x[1]*(L-14.4)
        print(" phi=%.1f p=%.2f: chi2 b=0: %.1f | b0 only: b0=%+.2f chi2 %.1f | b0+b1: b0=%+.2f b1=%+.2f chi2 %.1f | b(M2500=1.5e14)=%.2f b(5e14)=%.2f b(1e15)=%.2f | residual slope after fit %+.2f"%(phi,p,c0,r1.x[0],r1.fun,r2.x[0],r2.x[1],r2.fun,bf(np.log10(1.5e14)),bf(np.log10(5e14)),bf(15.0),np.polyfit(lM,np.log(fb/pr),1)[0]))
# what b1 is needed, with 1-sigma from chi2 profile (phi=0.7,p=0.5)
phi,p=.7,.5; best=minimize(lambda t: model(t,phi,p),[0.1,0.1],method="Nelder-Mead").x
grid=np.linspace(-0.1,0.6,29); prof=[]
for b1 in grid:
    r=minimize(lambda t: model([t[0],b1],phi,p),[best[0]],method="Nelder-Mead"); prof.append(r.fun)
prof=np.array(prof); ok=grid[prof<prof.min()+1]; print("\nprofile in b1 (phi=0.7,p=0.5): best b1=%.2f, 1-sigma [%.2f, %.2f]; b1=0 excluded by dchi2=%.1f"%(grid[np.argmin(prof)],ok.min(),ok.max(),model([minimize(lambda t: model([t[0],0],phi,p),[0.1],method='Nelder-Mead').x[0],0],phi,p)-prof.min()))
