#%%
import scf_pb
import numpy as np
# %%
scf_pb.D_eff_v(N=1000, sigma = 0.02, chi = np.linspace(0, 1,100), chi_PC = np.linspace(0, 1, 100), a0=0.18, a1=-0.09, particle_width=4, particle_height = 4, k_smooth = 4)
# %%
import inspect
# %%
scf_pb.D_v(N=1000, sigma = 0.02, chi = np.linspace(0, 1,100))
# %%
scf_pb.D_eff_v
# %%
