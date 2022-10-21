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
D_eff = scf_pb.D_eff_v(N=1000, sigma = 0.02, chi = np.linspace(0, 1, 10), chi_PC = np.linspace(0, 1, 50), a0=0.18, a1=-0.09, particle_width=4, particle_height = 4, k_smooth = 4, progressbar = True)
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
