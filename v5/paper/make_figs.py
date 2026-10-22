import io,contextlib,numpy as np,warnings; warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":150,"savefig.bbox":"tight"})
C={"fix":"#8c8c8c","v4":"#d95f02","v5":"#1b6ca8","obs":"#222222"}
OUT="paper/"
# ---- Fig 1: radial term
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('v4_vs_v5.py').read().split('print("X-COP 39')[0])
R_=lambda a: float(np.sqrt(np.mean(np.asarray(a)**2)))
m_fix=r_(K_(.5+A*xi,1.)); m_v4=r_((K_((.5+A*xi)*w**(-alpha),1.))); m_v5=r_(K_(.5+A*xi,dens**.5))
fig,ax=plt.subplots(1,3,figsize=(7.4,2.6),gridspec_kw={"width_ratios":[1.15,1,1]})
xs=np.arange(3); bw=.26
for k,(mm,lab,c) in enumerate(((m_fix,"n(ξ) at w=1",C["fix"]),(m_v4,"v4: w$^{-α(ξ)}$",C["v4"]),(m_v5,"this work: s=(Δ/Δ$_{500}$)$^{0.5}$",C["v5"]))):
    ax[0].bar(xs+(k-1)*bw,[R_(mm[:,j]) for j in range(3)],bw,color=c,label=lab)
ax[0].set_xticks(xs); ax[0].set_xticklabels(["0.5","1.0","1.5"]); ax[0].set_xlabel("r [Mpc]"); ax[0].set_ylabel("rms [dex]"); ax[0].set_title("(a) X-COP, 39 points",fontsize=9); ax[0].legend(frameon=False,fontsize=7)
for a,t,vals in ((ax[1],"(b) Vikhlinin: $f_{bar}(r_{2500})/f_{bar}(r_{500})$",[3.37,1.95,1.00,0.90]),(ax[2],"(c) X-COP: $f_{bar}(r_{200})/f_{bar}(r_{500})$",[0.42,0.97,0.97,1.03])):
    a.bar(range(4),vals,color=[C["fix"],C["v4"],C["v5"],C["obs"]]); a.set_xticks(range(4)); a.set_xticklabels(["no\nradial","v4","this\nwork","obs."],fontsize=7); a.axhline(1,color="k",lw=.5,ls=":"); a.set_title(t,fontsize=8)
    for i,v in enumerate(vals): a.text(i,v+.05,"%.2f"%v,ha="center",fontsize=7)
fig.tight_layout(); fig.savefig(OUT+"fig1_radial.pdf"); plt.close()
# ---- Fig 2: alpha scan
al=np.array([0,.25,.5,.75,1,1.25,1.5]); allen=[126.1,101.6,85.9,75.8,69.4,65.7,64.0]; mantz=[165.4,143.1,129.6,121.4,117.0,115.5,116.5]; gal=[.141,.113,.094,.097,.121]
fig,ax=plt.subplots(1,2,figsize=(6.6,2.6))
ax[0].plot(al,allen,"o-",color=C["v5"],label="Allen+08 (41 d.o.f.)"); ax[0].plot(al,mantz,"s-",color=C["v4"],label="Mantz+14 (39 d.o.f.)"); ax[0].set_xlabel(r"α in g$^\dagger\propto E(z)^{\alpha}$"); ax[0].set_ylabel("χ² (free normalisation)"); ax[0].legend(frameon=False,fontsize=7); ax[0].set_title("(a) clusters, r$_{2500}$",fontsize=9)
ax[1].plot([0,.25,.5,.75,1.0][:5],gal,"d-",color=C["obs"]); ax[1].set_xlabel("α"); ax[1].set_ylabel("rms of ln(V$_{pred}$/V$_{obs}$)"); ax[1].set_title("(b) KMOS3D medians (3 bins)",fontsize=9)
fig.tight_layout(); fig.savefig(OUT+"fig2_alpha.pdf"); plt.close()
# ---- Fig 3: Vikhlinin trend with mass
buf=io.StringIO()
with contextlib.redirect_stdout(buf): exec(open('vik_test.py').read())
fig,ax=plt.subplots(figsize=(3.4,2.7)); lmv=np.array(lm)
for k,c,mk,lab in (("p=0.50",C["v5"],"o","p=0.5"),("p(xi)=0.46(xi/0.4)^0.37",C["v4"],"s","p(ξ)=0.46(ξ/0.4)$^{0.37}$")):
    y=np.exp(np.array(res[k])); ax.scatter(lmv,y,c=c,marker=mk,s=18,label=lab)
    b=np.polyfit(lmv,np.log(y),1); xx=np.linspace(lmv.min(),lmv.max(),20); ax.plot(xx,np.exp(np.polyval(b,xx)),c=c,lw=.8)
ax.axhline(1,color="k",lw=.5,ls=":"); ax.set_xlabel(r"log M$_{2500}$ [M$_\odot$]"); ax.set_ylabel("$f_{bar}$ obs/pred at r$_{2500}$"); ax.legend(frameon=False,fontsize=7); fig.tight_layout(); fig.savefig(OUT+"fig3_mass_trend.pdf"); plt.close()
# ---- Fig 4: high-z BTFR
z=[.9,1.5,2.3]; Vo=[239,263,260]; V1=[215.0,290.0,302.9]; V0=[192.6,243.8,239.1]
fig,ax=plt.subplots(figsize=(3.4,2.7)); ax.plot(z,Vo,"ko-",label="KMOS3D median v$_{circ,max}$"); ax.plot(z,V1,"s--",color=C["v5"],label=r"RGD, g$^\dagger\propto H(z)$"); ax.plot(z,V0,"^--",color=C["fix"],label=r"RGD, g$^\dagger$ frozen")
ax.set_xlabel("z"); ax.set_ylabel("v [km/s]"); ax.legend(frameon=False,fontsize=7); fig.tight_layout(); fig.savefig(OUT+"fig4_highz.pdf"); plt.close()
print("ok", R_(m_fix),R_(m_v4),R_(m_v5), {k:np.polyfit(lmv,np.log(np.exp(np.array(res[k]))),1)[0] for k in res})
