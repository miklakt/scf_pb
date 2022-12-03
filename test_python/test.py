#%%
import scf_pb
import numpy as np
import time
# %%
#no threading 89.15948677062988 seconds
#1 worker 89.15948677062988 seconds
#4 workers 22.640379667282104 seconds
#10 workers 10.741898775100708 seconds
#16 workers 9.385826349258423 seconds
# access to cache 0.029059171676635742 seconds
start_time = time.time()
D_eff = scf_pb.D_eff_v(
    N=2000, sigma = 0.08, 
    chi = np.linspace(0, 1, 10), chi_PC = -1, 
    a0=0.18, a1=-0.09, 
    particle_width=np.geomspace(1,10), particle_height = "particle_width", 
    k_smooth = 1, 
    a="H", b = "H+1",
    progressbar = True
    )
print("--- %s seconds ---" % (time.time() - start_time))
# %%
scf_pb.phi(N=1000, chi=0.3, z=0, sigma = 0.02)
# %%
sigma = [0.01, 0.02, 0.03, 0.04]
# %%
scf_pb.phi_v(N=1000, chi=0.3, z=0, sigma = sigma)
# %%
chi = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
scf_pb.phi_v(N=1000, chi=chi, z=0, sigma = sigma)
# %%
scf_pb.mobility_factor_v(phi=0.1, d=10.0, k_smooth=1.0)
# %%
H = scf_pb.D(N=1000, sigma = 0.02, chi = 0.5)
z = np.linspace(0, H+100)
conc = scf_pb.conc_profile_fixed_source_v(
    N=1000, sigma = 0.02, 
    chi = 0.6, chi_PC = -0.5, 
    a0=0.18, a1=-0.09, 
    particle_width=8, particle_height = 8, 
    k_smooth = 1, 
    progressbar = True,
    source_dist=H+50,
    z=z
    )
plt.plot(z, conc)
# %%
start_time = time.time()
chi = np.linspace(0, 1, 50)
chi_PC = [0, -1.0,  -2.0, -4.0]
PC_sink = scf_pb.PC_perfect_sink_v(
    N=1000, sigma = 0.02, 
    chi_PC = chi_PC, chi = np.linspace(0, 1, 50),
    a0=0.18, a1=-0.09, 
    particle_width=8, particle_height = 8, 
    k_smooth = 1,
    l=0,
    progressbar = True
    )
PC_sink = scf_pb.PC_v(
    N=1000, sigma = 0.02, 
    chi_PC = chi_PC, chi = np.linspace(0, 1, 50),
    a0=0.18, a1=-0.09, 
    particle_width=8, particle_height = 8,
    progressbar = True
    )
print("--- %s seconds ---" % (time.time() - start_time))
#plt.plot(chi, PC)
plt.yscale("log")
[plt.plot(chi, PC_sink_, label = chi_PC_) for PC_sink_, chi_PC_ in zip(PC_sink, chi_PC)]
plt.legend()
# %%
