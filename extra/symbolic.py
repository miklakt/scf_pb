#%%
import sympy
Symbol = sympy.Symbol
# %%
phi = Symbol("\phi")
chi = Symbol("\chi")

Pi = -sympy.log(1-phi) - phi - phi**2*chi

# %%
phi_min_Pi = sympy.solve(Pi.diff(phi), phi)[1]
phi_min_Pi
# %%
Pi.subs(phi, phi_min_Pi).simplify()
# %%
