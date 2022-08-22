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
N=[1000.0]
sigma = [0.02]
chi_PC = np.linspace(0,-3)
chi_PS = np.linspace(0, 1)
a0 = [0.18]
a1 = [-0.09]
width = [8.0]
height = [8.0]
D_eff = scf_pb.D_eff_corrected_v(N, sigma, chi_PS, chi_PC, a0, a1, width, height)
#D_eff_corrected = scf_pb.D_eff_corrected_v(N, sigma, chi_PS, chi_PC, a0, a1, width, height)
D_eff_corrected=np.array(D_eff).reshape(len(chi_PS), len(chi_PC)).T
# %%
plt.imshow(np.log(D_eff_corrected), cmap='RdBu', extent=[min(chi_PS),max(chi_PS),min(chi_PC),max(chi_PC)])
plt.xlabel("$\chi_{PS}$")
plt.ylabel("$\chi_{PC}$")
cbar = plt.colorbar(extend='both')
cbar.set_label("$log(D_{eff}/D_s)$")
plt.clim(-5, 5)
#plt.contour(np.log(D_eff_corrected))
plt.title(f"d={width}")
plt.savefig(f"2/D_eff_corrected_im_{width[0]}_extended.pdf")
# %%
N=[1000.0]
sigma = [0.02]
chi_PC = [0, -0.5, -1.0, -1.5, -2.0]
chi_PS = np.linspace(0,1)
a0 = [0.18]
a1 = [-0.09]
width = [10.0]
height = [10.0]
D_eff = scf_pb.D_eff_v(N, sigma, chi_PS, chi_PC, a0, a1, width, height)
D_eff=np.array(D_eff).reshape(len(chi_PS), len(chi_PC)).T
D_eff_corrected = scf_pb.D_eff_corrected_v(N, sigma, chi_PS, chi_PC, a0, a1, width, height)
D_eff_corrected=np.array(D_eff_corrected).reshape(len(chi_PS), len(chi_PC)).T
# %%
fig, ax = plt.subplots()
for chi_pc_, d_eff, d_eff_corrected in zip(chi_PC, D_eff, D_eff_corrected):
    ax.plot(chi_PS, d_eff_corrected, label = chi_pc_)
    ax.plot(chi_PS, d_eff, color = LAST_USED_COLOR(), linewidth = 0.5, linestyle = ":")
    ax.set_xlabel("$\chi_{PS}$")
    ax.set_ylabel("$D_{eff}/D_s$")
ax.legend(title = "$\chi_{PC}$")
ax.set_yscale("log")
ax.set_title(f"d={width[0]}")
plt.ylim(1e-4, 1e4)
#plt.savefig(f"D_eff_corrected_{width[0]}_on_chi_PS.pdf")

# %%
N=[1000.0]
sigma = [0.02]
chi_PC = np.linspace(-2, 0)
chi_PS = [0, 0.25, 0.5, 0.75, 1.0]
a0 = [0.18]
a1 = [-0.09]
width = [4.0]
height = [4.0]
D_eff = scf_pb.D_eff_v(N, sigma, chi_PS, chi_PC, a0, a1, width, height)
D_eff=np.array(D_eff).reshape(len(chi_PS), len(chi_PC))
D_eff_corrected = scf_pb.D_eff_corrected_v(N, sigma, chi_PS, chi_PC, a0, a1, width, height)
D_eff_corrected=np.array(D_eff_corrected).reshape(len(chi_PS), len(chi_PC))
# %%
fig, ax = plt.subplots()
for chi_ps_, d_eff, d_eff_corrected in zip(chi_PS, D_eff, D_eff_corrected):
    ax.plot(chi_PC[::-1], d_eff_corrected, label = chi_ps_)
    ax.plot(chi_PC[::-1], d_eff, color = LAST_USED_COLOR(), linewidth = 0.5, linestyle = ":")
ax.set_xlabel("$\chi_{PC}$")
ax.set_ylabel("$D_{eff}/D_s$")
ax.legend(title = "$\chi_{PS}$")
ax.set_yscale("log")
ax.set_title(f"d={width[0]}")
plt.ylim(1e-4, 1e4)
plt.savefig(f"D_eff_corrected_{width[0]}_on_chi_PC.pdf")
plt.ylim(1e-1, 1e1)
plt.savefig(f"D_eff_corrected_{width[0]}_on_chi_PC_zoomed.pdf")
# %%
N=1000.0
sigma = 0.02
chi_PC = -2.8
chi_PS = [0, 0.25, 0.5, 0.75, 1.0]
a0 = 0.18
a1 = -0.09
d = 1/np.linspace(1/10, 1, 30)
D_eff_corrected = [[scf_pb.D_eff_corrected(N, sigma, chi_ps, chi_PC, a0, a1, d_, d_) for d_ in d] for chi_ps in chi_PS]
#%%
fig, ax = plt.subplots()
ax.plot(d, 1/d, color = "black")
for chi_ps_, d_eff in zip(chi_PS, D_eff_corrected):
    ax.plot(d, d_eff/d, label = chi_ps_)
ax.set_xlabel("$d$")
ax.set_ylabel(r"$D \frac{\eta_s}{k_B T}$")
ax.legend(title = "$\chi_{PS}$")
ax.set_title(f"{chi_PC=}")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_ylim(1e-2, 1e1)

#plt.savefig(f"D_eff_corrected_on_d_chi_PC_{chi_PC}.pdf")
# %%
fig, ax = plt.subplots()
ax.plot(1/(d*d), np.gradient(1/d,d), color = "black")
#for chi_ps_, d_eff in zip(chi_PS, D_eff_corrected):
#    ax.plot(1/(d*d), np.gradient(d_eff/d,d), label = chi_ps_)
ax.set_xlabel("$d^{-2}$")
ax.set_ylabel(r"$(D \frac{\eta_s}{k_B T})^\prime$")
ax.legend(title = "$\chi_{PS}$")
#ax.set_ylim(-1e1, 1e1)
#ax.set_title(f"{chi_PC=}")
#ax.set_xscale("log")
#ax.set_yscale("log")

#plt.savefig(f"D_eff_corrected_on_d_chi_PC_{chi_PC}.pdf")#
# %%
fig, ax = plt.subplots()
ax.plot(d, np.gradient(1/d,d), color = "black")
for chi_ps_, d_eff in zip(chi_PS, D_eff_corrected):
    ax.plot(d, np.gradient(d_eff/d,d), label = chi_ps_)
ax.set_xlabel("$d$")
ax.set_ylabel(r"$(D \frac{\eta_s}{k_B T})^\prime$")
ax.legend(title = "$\chi_{PS}$")
# %%
%matplotlib tk
from matplotlib.widgets import Slider, Button
fig, ax = plt.subplots()

ax.set_ylim(1e-2, 1e1)
ax.plot(d, 1/d, color = "black")

def f(chi_PC):
    D_eff_corrected = [[scf_pb.D_eff_corrected(N, sigma, chi_ps, chi_PC, a0, a1, d_, d_) for d_ in d] for chi_ps in chi_PS]
    return D_eff_corrected

lines = []
for chi_ps_, d_eff in zip(chi_PS, D_eff_corrected):
    line, = ax.plot(d, d_eff/d, label = chi_ps_)
    lines.append(line)

chi_PC_ax = plt.axes([0.2, 0.2, 0.6, 0.03])
chi_PC_slider = Slider(
    ax=chi_PC_ax,
    label="$\chi_{PC}$",
    valmin=-3,
    valmax=0,
    valinit=chi_PC,
)

def update(val):
    D_eff_corrected = f(chi_PC_slider.val)
    for chi_ps_, d_eff, line in zip(chi_PS, D_eff_corrected, lines):
        line.set_ydata(d_eff/d)
    fig.canvas.draw_idle()

chi_PC_slider.on_changed(update)
ax.legend(title = "$\chi_{PS}$")
ax.set_xlabel("$d$")
ax.set_ylabel(r"$D \frac{\eta_s}{k_B T}$")
ax.set_title(f"{chi_PC_slider.val=}")
ax.set_xscale("log")
ax.set_yscale("log")

# %%
