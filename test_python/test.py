#%%
import scf_pb
import numpy as np
# %%
scf_pb.D_eff_v(N=1000, sigma = 0.02, chi = np.linspace(0, 1,20), chi_PC = np.linspace(0, 1, 10), a0=0.18, a1=-0.09, particle_width=4, particle_height = 4, k_smooth = 4)
# %%
D = scf_pb.D(N=1000, sigma = 0.02, chi = 0)
phi = scf_pb.phi_v(N=1000, sigma = 0.02, chi = 0, z = np.arange(int(D+0.5)))
# %%
scf_pb.free_energy_external(phi=tuple(phi), chi=0.0, chi_PC=-1.0, a0=0.18, a1=-0.09, particle_width=4, particle_height=4, z=10)
# %%
scf_pb.D_eff(phi=tuple(phi), chi=0.0, chi_PC=-1.0, a0=0.18, a1=-0.09, particle_width=4, particle_height=4, a=2, b = D-2, k_smooth =4)
# %%
scf_pb.free_energy_all_v(N=1000, sigma = 0.02, chi=0.0, chi_PC=-1.0, a0=0.18, a1=-0.09, particle_width=4, particle_height=4, z=[10, 20, 30])
# %%
scf_pb._scf_pb.D_eff.__doc__
# %%

# %%
