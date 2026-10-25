"""PRE-REGISTERED (written before running).  Systematics-aware calibration of the radial parameter (v5: p; v6: q).
Points: every system at its measurement radius (ref52 D, 269) PLUS the r500 point of the paired samples (Vikhlinin 10, Ettori 9, groups 43; s=1 there, so they pin the sample offset).
Model: ln f_obs,ij = ln f_pred(theta) + delta_s + eps_ij,  delta_s ~ N(0, sd^2) per sample (Vikhlinin, Ettori, Allen, Groups-L(XMM), Groups-S(Chandra), XXL), eps ~ N(0, sf^2 + tau^2).
delta_s is marginalised analytically; tau profiled; theta on a grid.  Archived Theta (xi_c=.115, m=4), n(xi)=.5+.4068 xi, lam=.999, Kravtsov stars (phi as archived).
sd prior: 0.20 (central; the shared-cluster Chandra/XMM ratios 1.15-1.39 suggest 0.14-0.33), sensitivity 0.10, 0.30 and 'free' (sd=5: offsets unconstrained -> offset-free limit), and sd=0.01 (no offsets).
Outputs: best p and q with 1-sigma (dlnL=0.5), AIC (same number of free parameters: 1 + tau), LOSO log score (held-out sample's offset drawn from the prior), posterior offsets.
Decision rule for freezing (declared now): freeze v6 with q = best fit at sd=0.20 rounded to 2 decimals, and v5 with p likewise, as competing predictions; amplitude uncertainty sd=0.20 attached to every blind prediction."""
import io,contextlib,sys,numpy as np,warnings; warnings.filterwarnings("ignore")
sys.path.insert(0,'.')
with contextlib.redirect_stdout(io.StringIO()):
    import ref57 as Z
T=Z.T; D=Z.D; f500=Z.f500
XIC,MW=0.115,4
th=lambda xi: 1/(1+(XIC/xi)**MW)
pts=[]   # dict(samp, y, sf, rec, kind)
for d in D:
    gd=T.gd_of(d); xi=T.Gc*d['M5']/d['R5']**2/gd; dd=(d['M']/d['M5'])/(d['r']/d['R5'])**3; gr=(d['M']/d['M5'])/(d['r']/d['R5'])**2
    pts.append(dict(samp=d['samp'],y=np.log(d['fb']),sf=d['sf'],rec=d,xi=xi,dd=dd,gr=gr,r500=False))
    k=(d['samp'],d['name'])
    if k in f500:
        r5=dict(d,M=d['M5'],r=d['R5']); pts.append(dict(samp=d['samp'],y=np.log(f500[k]),sf=0.10,rec=r5,xi=xi,dd=1.0,gr=1.0,r500=True))
print("points:",len(pts)," r500 points:",sum(p['r500'] for p in pts))
samp=np.array([p['samp'] for p in pts]); y=np.array([p['y'] for p in pts]); sf=np.array([p['sf'] for p in pts])
xi=np.array([p['xi'] for p in pts]); dd=np.array([p['dd'] for p in pts]); gr=np.array([p['gr'] for p in pts]); T5=th(xi)
samples=sorted(set(samp)); idx={s:np.where(samp==s)[0] for s in samples}
def lnpred(law,par):
    out=np.empty(len(pts))
    for i,p in enumerate(pts):
        s=dd[i]**(par*T5[i]) if law=='v5' else gr[i]**(par*T5[i])
        f=T.fpred_s(p['rec'],s); out[i]=np.log(f) if np.isfinite(f) and f>0 else np.nan
    return out
G5=np.round(np.arange(0.30,1.001,0.025),3); G6=np.round(np.arange(0.60,1.601,0.025),3)
RES={('v5',g):y-lnpred('v5',g) for g in G5}; RES.update({('v6',g):y-lnpred('v6',g) for g in G6})
def ll_sample(r,s_,tau,sd):
    v=s_**2+tau**2; w=1/v; A=w.sum(); B=(w*r).sum(); k=sd**2
    return -0.5*(np.sum(np.log(2*np.pi*v))+np.sum(w*r**2)-k*B**2/(1+k*A)+np.log(1+k*A))
def ll_all(r,tau,sd,samps=None):
    tot=0
    for s in (samples if samps is None else samps):
        ii=idx[s]; rr=r[ii]; 
        if np.any(~np.isfinite(rr)): return -1e9
        tot+=ll_sample(rr,sf[ii],tau,sd)
    return tot
from scipy.optimize import minimize_scalar
def prof(law,g,sd,samps=None):
    r=RES[(law,g)]; o=minimize_scalar(lambda lt:-ll_all(r,np.exp(lt),sd,samps),bounds=(-6,1),method='bounded'); return -o.fun,np.exp(o.x)
def scan(law,sd):
    G=G5 if law=='v5' else G6; L=np.array([prof(law,g,sd)[0] for g in G]); j=int(np.argmax(L))
    # parabola refine
    if 0<j<len(G)-1:
        c=np.polyfit(G[j-1:j+2],L[j-1:j+2],2); gb=-c[1]/(2*c[0]); Lb=np.polyval(c,gb); se=np.sqrt(-1/(2*c[0]))
    else: gb,Lb,se=G[j],L[j],np.nan
    return gb,se,Lb,L,G
print("\nsd_delta   | v5: p_best ± se   lnL   | v6: q_best ± se   lnL   | dAIC(v6-v5) | lnL gain")
store={}
for sd in (0.01,0.10,0.20,0.30,5.0):
    a=scan('v5',sd); b=scan('v6',sd); store[sd]=(a,b)
    print("%-10s | %.3f ± %.3f  %8.1f | %.3f ± %.3f  %8.1f | %+8.1f | %+.1f"%("free" if sd==5 else sd,a[0],a[1],a[2],b[0],b[1],b[2],-2*(b[2]-a[2]),b[2]-a[2]))
# fixed archived vs best, at sd=0.2, and LOSO
sd=0.20
def lnl_at(law,par):
    G=G5 if law=='v5' else G6; g=G[np.argmin(abs(G-par))]; return prof(law,g,sd)[0]
print("\nat sd=0.20: lnL(v5,p=0.50)=%.1f  lnL(v6,q=1.00)=%.1f  | best: v5 %.1f  v6 %.1f"%(lnl_at('v5',.5),lnl_at('v6',1.0),store[sd][0][2],store[sd][1][2]))
print("\nLOSO (theta & tau fit on the other samples, held-out offset from prior sd=0.20): held-out log score")
print("%-18s  v5(p_chosen)   v6(q_chosen)   v5(p=.5)  v6(q=1)   chosen p,q"%"sample")
tot={'v5':0,'v6':0,'v5f':0,'v6f':0}
for s in samples:
    others=[x for x in samples if x!=s]; row={}
    for law,G in (('v5',G5),('v6',G6)):
        L=[prof(law,g,sd,others)[0] for g in G]; gb=G[int(np.argmax(L))]; tau=prof(law,gb,sd,others)[1]
        row[law]=(gb,ll_all(RES[(law,gb)],tau,sd,[s]))
        gf=0.5 if law=='v5' else 1.0; gi=G[np.argmin(abs(G-gf))]; tf=prof(law,gi,sd,others)[1]; row[law+'f']=ll_all(RES[(law,gi)],tf,sd,[s])
    for k in ('v5','v6'): tot[k]+=row[k][1]
    tot['v5f']+=row['v5f']; tot['v6f']+=row['v6f']
    print("%-18s  %8.1f       %8.1f      %8.1f  %8.1f   %.2f, %.2f"%(s,row['v5'][1],row['v6'][1],row['v5f'],row['v6f'],row['v5'][0],row['v6'][0]))
print("%-18s  %8.1f       %8.1f      %8.1f  %8.1f"%("TOTAL",tot['v5'],tot['v6'],tot['v5f'],tot['v6f']))
# posterior offsets at best theta, sd=0.2
print("\nposterior offsets delta_s (ln f) at sd=0.20, best theta (positive = observed above prediction):")
for law,(gb,se,Lb,L,G) in zip(('v5','v6'),store[0.20]):
    g=G[np.argmin(abs(G-gb))]; r=RES[(law,g)]; tau=prof(law,g,0.20)[1]; out=[]
    for s in samples:
        ii=idx[s]; v=sf[ii]**2+tau**2; w=1/v; A=w.sum(); B=(w*r[ii]).sum(); k=0.20**2; out.append("%s %+.2f"%(s.split(' ')[0].split('(')[0],k*B/(1+k*A)))
    print("  %s (par %.3f, tau %.3f): "%(law,g,tau)+"; ".join(out))
