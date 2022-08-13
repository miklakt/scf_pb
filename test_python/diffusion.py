#%%
import numpy as np
import matplotlib.pyplot as plt
#matplotlib settings
LAST_USED_COLOR = lambda: plt.gca().lines[-1].get_color()
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
# %%
import scf_pb
# %%
chi = 0.5
N=1000.0
sigma = 0.02
R = 300.0
D = scf_pb.D(N, sigma, chi)
# %%
z = np.linspace(0, D)
phi =[scf_pb.phi(N, sigma, chi, z_) for z_ in z]
# %%
plt.plot(z, phi)
# %%
CHI_PS = np.linspace(0,1)
D = [scf_pb.D(N, sigma, chi) for chi in CHI_PS]
# %%
plt.plot(CHI_PS, D)
# %%
def diffusion_factor(phi, d, n=1):
    if phi == 0:
        return 1
    eps = 1/phi
    a = eps**2/d**2
    return a/((1+a**(n)))**(1/n)

def diffusion_piecewise(phi, d):
    eps = 1/phi
    if phi == 0:
        return 1
    if d<eps:
        return 1
    else:
        return eps**2/d**2
# %%
phi = np.linspace(0,1)
d=4
factor = [diffusion_piecewise(phi_, d) for phi_ in phi]
factor_smooth = [diffusion_factor(phi_, d, n=2) for phi_ in phi]
# %%
fig, ax = plt.subplots(nrows=2, sharey=True)
ax[0].plot(phi, factor, label = "piecewise")
ax[0].plot(phi, factor_smooth, label = r"$\frac{x}{{(1-x^k)}^{1/k}}$")
#ax[0].set_yscale("log")
ax[0].set_xlabel("$\phi$")
ax[0].set_ylabel("$D/D_s$")
ax[0].legend()

ax[1].plot(phi**(-1),factor)
ax[1].plot(phi**(-1),factor_smooth)
ax[1].set_xlabel("$\phi^{-1}$")
ax[1].set_ylabel("$D/D_s$")
ax[1].set_xscale("log")

plt.tight_layout()
# %%
chi_PC = -1.5
CHI_PS = np.linspace(0,1)
a0 = 0.18
a1 = -0.09
width = 4.0
height = 4.0
D_eff = [scf_pb.D_eff(N, sigma, chi, chi_PC, a0, a1, width, height) for chi in CHI_PS]
D_eff_corrected = [scf_pb.D_eff_corrected(N, sigma, chi, chi_PC, a0, a1, width, height) for chi in CHI_PS]
# %%
fig, ax = plt.subplots()
ax.plot(CHI_PS, D_eff)
ax.plot(CHI_PS, D_eff_corrected, label = "corrected")
ax.set_xlabel("$\chi_{PS}$")
ax.set_ylabel("$D_{eff}/D_s$")
ax.legend()
ax.set_yscale("log")

# %%
