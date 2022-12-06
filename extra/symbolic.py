#%%
import sympy
Symbol = sympy.Symbol
# %%
phi = Symbol("\phi")
chi = Symbol("\chi")
chi_PC = Symbol("\chi_{PC}")
chi_crit =  Symbol("\chi_{crit}")
d = Symbol("d")
a0 = Symbol("a_0")
a1 = Symbol("a_1")
Pi = -sympy.log(1-phi) - phi - phi**2*chi
gamma = (chi_PC - chi*(1-phi) - chi_crit)*(a0*phi+a1*phi**2)
# %%
V = sympy.pi*d**3/6
S = sympy.pi*d**2
# %%
F = Pi*V+gamma*S
# %%
dF_dd = sympy.diff(F, d)
# %%
d_crit=sympy.solve(dF_dd/d,d)[0]
# %%
chi_PC_min=sympy.solve(d_crit, chi_PC)[0]
# %%
from sympy.plotting import plot
x = Symbol("x")
plot(d_crit
    .subs(chi_PC, -1)
    .subs(phi, 0.4)
    .subs(chi_crit, -1.1)
    .subs(a0, 0.18)
    .subs(a1, -0.09)
    ,(chi, 0,1))
# %%
