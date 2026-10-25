"""PRE-REGISTERED follow-up (declared after seeing ref52 results; exploratory, so uses AIC+LOSO only to *rank*, not to claim):
ref52 found the best one-parameter ratio is gr=g(r)/g500 with q~1.1 (better than dd=Delta/Delta500 by dAIC~55). Questions:
 (1) is q=1 (parameter-free, s=g(r)/g500) acceptable?  (2) two-parameter family s=dd^a*rr^b (gr: a=b=1; dd: b=0).  (3) does an XMM offset or an offset c remain?  (4) which samples drive it (per-sample LOSO).
With s=g/g500 the kernel argument becomes f*xi (independent of r): reported as an implication, not assumed."""
import sys,numpy as np; sys.path.insert(0,'.')
import ref52 as T
from ref52 import *
def lns2(model,th,idx):
    if model=='gr1': return np.log(W['gr'][idx])
    if model=='grc': return th[0]*np.log(W['gr'][idx])+th[1]
    if model=='grx': return th[0]*np.log(W['gr'][idx])+th[1]*xmm[idx]
    if model=='dr': return th[0]*np.log(W['dd'][idx])+th[1]*np.log(W['rr'][idx])
    if model=='drx': return th[0]*np.log(W['dd'][idx])+th[1]*np.log(W['rr'][idx])+th[2]*xmm[idx]
    return None
old=T.lns_fn
def new(model,th,idx):
    r=lns2(model,th,idx); return old(model,th,idx) if r is None else r
T.lns_fn=new; lns_fn=new
T.npar.update({'gr1':0,'grc':2,'grx':2,'dr':2,'drx':3})
# resid/nll/fit in ref52 reference module globals lns_fn -> patched via T.lns_fn
base=T.nll('q_gr',[1.125],np.log(0.234),allidx)
models=[('q_gr',[1.1]),('gr1',[]),('grc',[1.,0.]),('grx',[1.,0.]),('q_dd',[.5]),('dr',[.5,.5]),('drx',[.5,.5,0.])]
res={}
print("model   k  tau    nll    AIC    LOSO   params")
for m,th0 in models:
    th,tau,v=T.fit(m,allidx,th0); k=T.npar.get(m,1); aic=2*v+2*(k+1)
    lo=0; per={}
    for sm in samples:
        te=np.where(samp==sm)[0]; tr=np.where(samp!=sm)[0]
        tht,taut,_=T.fit(m,tr,th0); per[sm]=-T.nll(m,tht,np.log(taut),te); lo+=per[sm]
    res[m]=per
    print("%-7s %d  %.3f  %6.1f %6.1f %7.1f  %s"%(m,k,tau,v,aic,lo,np.round(th,3)))
print("\nper-sample LOSO log score: sample | "+" ".join("%8s"%m for m,_ in models))
for sm in samples: print("%-18s| "%sm+" ".join("%8.1f"%res[m][sm] for m,_ in models))
