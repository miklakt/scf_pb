#%%
import matplotlib.pyplot as plt
import numpy as np
import scf_pb

# %%
r=5
d = r*2
a0 = 0.18
a1 = -0.09
chi = 1
sigma = 0.08
N = 2000

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
chi_PC=-3
chi_ads = chi_PC - chi*(1-phi)
chi_crit = 6*np.log(5/6)
gamma = (chi_ads-chi_crit)*(a0*phi+a1*phi**2)
fe_approx =  4/3*np.pi*r**3*Pi + 4*np.pi*r**2*gamma
#fe_approx =  np.pi*d**3/4*Pi# +3/2*np.pi*d**2*gamma
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
d=np.linspace(1,10)
fe = scf_pb.free_energy_v(
    N=N, sigma=sigma, chi=chi,
    chi_PC = chi_PC,
    a0=a0, a1=a1, 
    particle_width = d, particle_height = "particle_width", 
    z = "H"
    )
plt.plot(d, fe)
# %%
d=np.linspace(1,10)
D_eff = scf_pb.D_eff_v(
    N=N, sigma=sigma, chi=chi,
    chi_PC = chi_PC,
    a0=a0, a1=a1, 
    particle_width = d, particle_height = "particle_width", 
    k_smooth = 1,
    a=0, b=["H"],
    progressbar=True
    )
plt.plot(d, D_eff)
# %%
d=np.linspace(1,10)
D_eff = scf_pb.D_eff_v(
    N=N, sigma=sigma, chi=chi,
    chi_PC = chi_PC,
    a0=a0, a1=a1, 
    particle_width = d, particle_height = "particle_width", 
    k_smooth = 1,
    a="H", b="H+10",
    progressbar=True
    )
plt.plot(d, D_eff)
plt.xscale("log")
plt.yscale("log")
# %%
