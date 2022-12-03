#%%
import numpy as np
import matplotlib.pyplot as plt
#matplotlib settings
LAST_USED_COLOR = lambda: plt.gca().lines[-1].get_color()
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
#%%
import scf_pb
# %%
def Pi_py(phi, chi):
    try:
        Pi = -np.log(1-phi) - phi - chi*phi**2
    except:
        Pi=0
    return Pi
#%%
N=2000
sigma = 0.08
chi = np.arange(0, 1.1, 0.1)
D = scf_pb.D_v(N=N, sigma=sigma, chi=chi)+100
z = [np.linspace(0, D_) for D_ in D]
phi = [scf_pb.phi_v(N=N, sigma=sigma, chi=chi_, z = z_) for z_, chi_ in zip(z, chi)]
Pi = [scf_pb.Pi_v(N=N, sigma=sigma, chi=chi_, z = z_) for z_, chi_ in zip(z, chi)]
Pi2 = [[Pi_py(phi__, chi_) for phi__ in phi_] for phi_, chi_ in zip(phi, chi)]
# %%
[plt.plot(z_, phi_, label = f"{chi_:.1f}") for z_, phi_, chi_ in zip(z, phi, chi)]
plt.legend(title= "$\chi$")
plt.xlabel("z")
plt.ylabel("$\phi$")
# %%
[plt.plot(z_, Pi_, label = f"{chi_:.1f}") for z_, Pi_, chi_ in zip(z, Pi, chi)]
plt.legend(title= "$\chi$")
plt.xlabel("z")
plt.ylabel("$\phi$")
#%%
phi_D =scf_pb.phi_v(N=N, sigma=sigma, chi=chi, z = "H")
plt.plot(chi, phi_D)
plt.legend(title= "$\chi$")
plt.xlabel("\chi")
plt.ylabel("$\phi_D$")
#%%
phi = np.linspace(0,0.5)
Pi=np.vectorize(Pi_py)(phi, 0.7)
plt.plot(phi, Pi)
plt.axhline(0)
# %%
import scipy.optimize
# %%
scipy.optimize.fsolve(lambda x: Pi_py(x, 0.7), 0.28)
# %%
