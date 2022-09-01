#%%
import scf_pb
# %%
scf_pb.D_eff(N=1000, sigma = 0.02, chi = [0, 0.1], chi_PC = [-2, -2.5], a0=0.18, a1=-0.09, particle_width=4, particle_height = 4, k_smooth = 4)
# %%
import inspect
inspect.signature(scf_pb.D_eff_cxx)
# %%
doc =scf_pb._scf_pb.D_eff.__doc__
# %%
sig = doc.split("\n")[0]
# %%
import inspect
inspect.signature(scf_pb.D_eff_cxx)
# %%
scf_pb.D_eff_cxx.__doc__
#%%
scf_pb.D_eff
# %%
