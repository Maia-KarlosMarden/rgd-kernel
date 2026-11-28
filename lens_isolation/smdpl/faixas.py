import json, sys, numpy as np
d=json.load(open(sys.argv[1]))
bins=['logMb_9.0-10.5','logMb_10.5-10.8','logMb_10.8-11.1','logMb_11.1-11.5']
def comb(a):
    w=np.array([1/e**2 for _,e in a]); v=np.array([x for x,_ in a])
    return (w*v).sum()/w.sum(), w.sum()**-0.5
for tag,r in d.items():
    s=r['amostras']; b=[(s[n]['eps_0.3_1Mpc'],s[n]['erro']) for n in bins]
    lo,el=comb(b[:2]); hi,eh=comb(b[2:])
    print(f"{tag:7s}", ' '.join(f"{v:.2f}±{e:.2f}" for v,e in b), f"| baixa-alta={lo-hi:+.2f}±{np.hypot(el,eh):.2f}",
          "| isoladas:", [s[n]['n_isoladas'] for n in bins])
