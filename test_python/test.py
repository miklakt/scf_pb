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
scf_pb.D_eff_v(N=1000, sigma = 0.02, chi = np.linspace(0, 1, 50), chi_PC = np.linspace(0, 1, 50), a0=0.18, a1=-0.09, particle_width=4, particle_height = 4, k_smooth = 4, progressbar = True, max_workers = 4)
print("--- %s seconds ---" % (time.time() - start_time))