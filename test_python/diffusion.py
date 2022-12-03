#%%
import matplotlib.pyplot as plt
import numpy as np
import scf_pb

# %%
d = 10
a0 = 0.18
a1 = -0.09
chi = 0.6
chi_PC =-1
sigma = 0.02
N = 1000

H = scf_pb.D(N=N, sigma= sigma, chi = chi)
z = np.linspace(0, H*1.5)
#%%
phi = scf_pb.phi_v(N=N, sigma=sigma, chi = chi, z=z).astype(float)
plt.plot(z, phi)
# %%
Pi = scf_pb.Pi_v(N=N, sigma=sigma, chi = chi, z=z)
Pi2 = -np.log(1-phi) - phi - chi*phi**2
plt.plot(z, Pi)
plt.plot(z, Pi2)
#%%
def volume(w, h = None):
    if h is None:
        h=w
    return np.pi*w**2/4*h

def surface(w, h = None):
    if h is None:
        h=w
    return np.pi*w*h+np.pi*w**2/2

def gamma(a1, a2, chi_PS, chi_PC, phi):
    chi_crit = 6*np.log(5/6)
    chi_ads = chi_PC - chi_PS*(1-phi)
    gamma = (chi_ads - chi_crit)*(a1*phi+a2*phi**2)
    return gamma

def Pi(phi, chi_PS):
    return -np.log(1-phi) - phi - chi_PS*phi**2

def surface_free_energy(phi, a1, a2, chi_PS, chi_PC, w, h=None):
    return surface(w, h)*gamma(a1, a2, chi_PS, chi_PC, phi)

def volume_free_energy(phi, chi_PS, w, h=None):
    return Pi(phi, chi_PS)*volume(w,h)

def free_energy_penalty_phi(phi, a1, a2, chi_PS, chi_PC, w, h=None):
    return surface_free_energy(phi, a1, a2, chi_PS, chi_PC, w, h)+volume_free_energy(phi, chi_PS, w, h)
#%%
fe_approx = [free_energy_penalty_phi(phi_, a0, a1, chi, chi_PC, d) for phi_ in phi]
fe = scf_pb.free_energy_v(
    N=N, sigma=sigma, chi=chi,
    chi_PC = chi_PC,
    a0=a0, a1=a1, 
    particle_width = d, particle_height = "particle_width", 
    z = z
    )
plt.plot(z, fe, label ="fe")
plt.plot(z, fe_approx,  label ="approx")
plt.legend()
# %%
