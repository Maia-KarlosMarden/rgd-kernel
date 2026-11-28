"""Teste do argumento contra viscosidade negativa como fonte de gravidade extra. Sem dados, sem ajuste."""
import numpy as np, sympy as sp

t,r,th,ph=X=sp.symbols('t r theta phi')
Phi,Lam=sp.Function('Phi')(r),sp.Function('Lambda')(r)
g=sp.diag(-sp.exp(2*Phi),sp.exp(2*Lam),r**2,r**2*sp.sin(th)**2); gi=g.inv()
Gam=[[[sum(gi[a,d]*(sp.diff(g[d,b],X[c])+sp.diff(g[d,c],X[b])-sp.diff(g[b,c],X[d])) for d in range(4))/2
       for c in range(4)] for b in range(4)] for a in range(4)]

print("== 1. Meio em repouso num espaco-tempo estatico (u paralelo ao vetor de Killing)")
uu=sp.Matrix([sp.exp(-Phi),0,0,0]); ud=g*uu
D=sp.Matrix(4,4,lambda a,b: sp.diff(ud[b],X[a])-sum(Gam[c][a][b]*ud[c] for c in range(4)))  # D[a,b]=nabla_a u_b
acc=sp.Matrix([sum(uu[a]*D[a,b] for a in range(4)) for b in range(4)])
h=g+ud*ud.T; theta=sp.simplify(sum(gi[a,b]*D[a,b] for a in range(4) for b in range(4)))
sig=sp.Matrix(4,4,lambda a,b: sp.simplify((D[a,b]+D[b,a])/2+(acc[b]*ud[a]+acc[a]*ud[b])/2-theta*h[a,b]/3))
print("theta =",theta,"| sigma nulo:",sig==sp.zeros(4),"| aceleracao_r =",sp.simplify(acc[1]))

print("\n== 2. Brecha: meio em fluxo radial v=-sqrt(2GM/r) (limite newtoniano)")
GM=sp.symbols('GM',positive=True); rr=sp.symbols('r',positive=True); v=-sp.sqrt(2*GM/rr)
print("theta =",sp.simplify(sp.diff(v,rr)+2*v/rr),"| sigma_rr =",sp.simplify(sp.Rational(2,3)*(sp.diff(v,rr)-v/rr)))

print("\n== 3. O que cada tensao gera (Einstein linearizado)")
x,y,z,e=sp.symbols('x y z e'); Y=[t,x,y,z]
F=sp.Function('F')(x,y,z); P=sp.Function('P')(x,y,z)
gl=sp.diag(-(1+2*e*F),1-2*e*P,1-2*e*P,1-2*e*P)
eta=sp.diag(-1,1,1,1)
lin=lambda ex: sp.diff(ex,e).subs(e,0)
G1=[[[sum(eta[a,d]*lin(sp.diff(gl[d,b],Y[c])+sp.diff(gl[d,c],Y[b])-sp.diff(gl[b,c],Y[d])) for d in range(4))/2
      for c in range(4)] for b in range(4)] for a in range(4)]
Ric=lambda b,c: sum(sp.diff(G1[a][b][c],Y[a])-sp.diff(G1[a][b][a],Y[c]) for a in range(4))
R=sum(eta[a,a]*Ric(a,a) for a in range(4))
Ein=lambda a: sp.simplify(Ric(a,a)-eta[a,a]*R/2)
lap=lambda f: sum(sp.diff(f,q,2) for q in (x,y,z))
print("G00 - 2 lap(Psi) =",sp.simplify(Ein(0)-2*lap(P)))
print("traco espacial - 2 lap(Phi-Psi) =",sp.simplify(sum(Ein(i) for i in (1,2,3))-2*lap(F-P)))
print("=> lap Psi = 4piG rho ; lap Phi = 4piG(rho + 3p/c2). Pressao p sem energia: dinamica recebe 3p/c2, lente (Phi+Psi)/2 recebe metade.")

print("\n== 4. Ordens de grandeza (exemplo redondo, nao e dado: massa extra 1e15 Msun em 1 Mpc)")
G=6.674e-11; c=2.998e8; Mpc=3.0857e22; Msun=1.989e30; H=67.4e3/Mpc
M=1e15*Msun; Rc=Mpc; rho_x=M/(4/3*np.pi*Rc**3); p=rho_x*c**2/3
rho_c=3*H**2/(8*np.pi*G); zeta=p/(3*H)
print(f"rho_extra={rho_x:.2e} kg/m3  p necessaria={p:.2e} Pa  |zeta| (theta=3H)={zeta:.2e} Pa s")
print(f"pressao homogenea 3H|zeta| / (rho_crit c2) = {3*H*zeta/(rho_c*c**2):.0f}")
gN=G*M/Rc**2
print(f"gradiente de pressao / peso do meio ~ c2/(g r) = {c**2/(gN*Rc):.1e}")

print("\n== 5. Estabilidade, modo transversal: rho s(1+tau s) + eta k2 = 0, com eta<0 (unidades arbitrarias)")
for tau in (0.0,1.0):
  for k in (1,10,100,1000):
    s=k**2 if tau==0 else max(np.roots([tau,1,-k**2]).real)
    print(f"tau={tau:g} k={k:5d}  taxa de crescimento={s:.4g}")
k_,nu,kap=sp.symbols('k nu kappa',positive=True)
print("com estabilizador -kappa k^4, pico em k =",sp.solve(sp.diff(nu*k_**2-kap*k_**4,k_),k_))

print("\n== 6. Ondas gravitacionais: amplitude x exp(-8 pi G eta D/c^3), D=40 Mpc")
D40=40*Mpc
for frac in (0.1,0.2):
    print(f"desvio {frac:.0%}: |eta| < {np.log(1+frac)*c**3/(8*np.pi*G*D40):.2e} Pa s")
