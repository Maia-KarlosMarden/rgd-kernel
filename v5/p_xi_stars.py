import io,contextlib,numpy as np
from scipy.optimize import minimize
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('hse_bias.py').read().split('print("Allen+2008')[0])
phi=.7
stG=np.array([phi*Mstar_tot(r['M5']/Msun)*Msun/(r['f25']*r['M25']) for r in grows]); xiG=np.array([r['xi'] for r in grows])
stA=phi*Mstar_tot(M5/Msun)*Msun/(fo*M25)
def pfun(x,th): return th[0]*(x/0.4)**th[1]
def cost(th,use):
    if th[0]<=0.05: return 1e9
    t=0
    if "G" in use: rr=np.array([gresid(r,stG[i],1,p=pfun(xiG[i],th)) for i,r in enumerate(grows)]); t+=len(rr)*np.log(np.mean(rr**2))
    if "A" in use:
        pr=np.array([fsolve(m,r,x,c_*Hz(zz)/(2*np.pi),pfun(x,th)) for m,r,x,zz in zip(M25,r25,xi,z)]); fb=fo*(1+stA); d=np.log(fb/pr); t+=len(d)*np.log(np.mean(d**2))
    return t
for use,nm in ((("G","A"),"groups+Allen"),(("A",),"Allen only"),(("G",),"groups only")):
    best=None
    for s in ([.5,0.],[.5,.5],[.5,-.5]):
        r=minimize(cost,s,args=(use,),method="Nelder-Mead",options=dict(xatol=1e-3,fatol=1e-3))
        if best is None or r.fun<best.fun: best=r
    c0=min(minimize(cost,[.5,0.],args=(use,),method="Nelder-Mead").fun,1e9); 
    rc_=minimize(lambda t: cost([t[0],0.0],use),[.5],method="Nelder-Mead")
    th=best.x
    pr=np.array([fsolve(m,r,x,c_*Hz(zz)/(2*np.pi),pfun(x,th)) for m,r,x,zz in zip(M25,r25,xi,z)]); fb=fo*(1+stA)
    print(" fit on %-13s: p0=%.2f q=%+.2f (p at xi=0.15: %.2f, 0.25: %.2f, 0.4: %.2f, 0.6: %.2f, 0.8: %.2f) | -2lnL %.1f vs constant-p %.1f (p=%.2f) | Allen mass slope: obs %+.2f model %+.2f ; Allen obs/pred %.2f"%(nm,th[0],th[1],*[pfun(x,th) for x in (.15,.25,.4,.6,.8)],best.fun,rc_.fun,rc_.x[0],np.polyfit(lM,np.log(fb),1)[0],np.polyfit(lM,np.log(pr),1)[0],np.median(fb/pr)))
