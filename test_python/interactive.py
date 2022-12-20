#%%
import scf_pb
import matplotlib.pyplot as plt
import numpy as np

import matplotlib
#matplotlib.use('TkAgg')
#%%
N =2000
sigma = 0.08
a0 = 0.18
a1=-0.09
def get_phi_profile(N, sigma, chi):
    D = scf_pb.D(N=N, sigma = sigma, chi = chi)+10
    z_arr = np.linspace(0,D, 500)
    phi = scf_pb.phi_v(N=N, sigma=sigma, chi=chi, z = z_arr)
    z_arr = np.append(z_arr, z_arr[-1])
    phi = np.append(phi, 0)
    return z_arr, phi

def get_free_energy_profile(N, sigma, chi, chi_PC, d):
    D = scf_pb.D(N=N, sigma = sigma, chi = chi)+10
    z_arr = np.linspace(0,D, 500)
    fe_all = np.array(scf_pb.free_energy_all_v(N=N, sigma = sigma, chi = chi, chi_PC = chi_PC, a0=a0, a1=a1, particle_width = d, particle_height = d, z = z_arr).tolist())
    fe, fe_osm, fe_surf = fe_all.swapaxes(-1, 0)
    #fe_osm = scf_pb.free_energy_osm_v(N=N, sigma = sigma, chi = chi, particle_width = d, particle_height = d, z = z_arr)
    #fe_surf = scf_pb.free_energy_surf_v(N=N, sigma = sigma, chi = chi, chi_PC = chi_PC, a0=a0, a1=a1, particle_width = d, particle_height = d, z = z_arr)
    z_arr = np.append(z_arr, z_arr[-1])
    fe = np.append(fe, 0)
    fe_surf = np.append(fe_surf, 0)
    fe_osm = np.append(fe_osm, 0)
    return z_arr, fe, fe_osm, fe_surf
# %%
%matplotlib ipympl
from matplotlib.widgets import Slider, Button
fig, ax = plt.subplots(nrows=3, sharex=True)
plt.subplots_adjust(bottom=0.3)

chi_slider = Slider(
    ax=plt.axes([ax[-1].get_position().x0,
                     ax[-1].get_position().y0-0.1,
                     ax[-1].get_position().width, 0.03]),
    label="$\chi_{PS}$",
    valmin=0,
    valmax=1,
    valinit=0.5,
)

chi_PC_slider = Slider(
    ax=plt.axes([ax[-1].get_position().x0,
                     ax[-1].get_position().y0-0.125,
                     ax[-1].get_position().width, 0.03]),
    label="$\chi_{PC}$",
    valmin=-3,
    valmax=-2,
    valinit=0,
)

d_slider = Slider(
    ax=plt.axes([ax[-1].get_position().x0,
                     ax[-1].get_position().y0-0.15,
                     ax[-1].get_position().width, 0.03]),
    label="$d$",
    valmin=1,
    valmax=16,
    valinit=8,
)

z, phi = get_phi_profile(N, sigma, chi_slider.val)
#z = np.arange(0,200, 0.5)
#phi = scf_pb.phi_v(N=N, sigma=sigma, chi=chi_PC_slider.val, z=z)
_, fe, fe_osm, fe_surf = get_free_energy_profile(N, sigma, chi_slider.val, chi_PC_slider.val, d_slider.val)
mobility_factor = scf_pb.D_z_v(N=N, sigma=sigma, chi = chi_slider.val, chi_PC= chi_PC_slider.val, particle_width = d_slider.val, particle_height =  d_slider.val, k_smooth = 1, z=z, a0=a0, a1=a1)
conc_profile = scf_pb.conc_profile_v(N=N, sigma=sigma, chi = chi_slider.val, chi_PC= chi_PC_slider.val, particle_width = d_slider.val, particle_height =  d_slider.val, k_smooth = 1, z=z, a0=a0, a1=a1, l=10)

line, = ax[0].plot(z, phi, label = "$\phi$")
line2, = ax[1].plot(z, fe, label = "$F_{tot}$")
line21, = ax[1].plot(z, fe_osm, label = "$F_{osm}$")
line22, = ax[1].plot(z, fe_surf, label = "$F_{surf}$")
line31, = ax[2].plot(z, conc_profile, label = "$F_{surf}$")

twinax = ax[0].twinx()
twinax.set_ylim(0,1.1)
line3, = twinax.plot(z, mobility_factor[::-1], color = "red", label = "mobility factor")
ax[0].plot([],[], color = "red", label = "mobility factor")

phi_max = 0
z_max = 0

ax[2].set_ylim(0,np.max(conc_profile))

def update(val):
    global phi_max
    global z_max
    chi_PS = chi_slider.val

    z, phi = get_phi_profile(N, sigma, chi_slider.val)
    _, fe, fe_osm, fe_surf = get_free_energy_profile(N, sigma, chi_slider.val, chi_PC_slider.val, d_slider.val)
    mobility_factor = scf_pb.D_z_v(N=N, sigma=sigma, chi = chi_slider.val, chi_PC= chi_PC_slider.val, particle_width = d_slider.val, particle_height =  d_slider.val, k_smooth = 1, z=z, a0=a0, a1=a1)
    conc_profile = scf_pb.conc_profile_v(N=N, sigma=sigma, chi = chi_slider.val, chi_PC= chi_PC_slider.val, particle_width = d_slider.val, particle_height =  d_slider.val, k_smooth = 1, z=z, a0=a0, a1=a1, l=10)
    phi_max = max(phi_max, max(phi))


    z_max = max(z_max, max(z))

    line.set_xdata(z)
    line.set_ydata(phi)

    line2.set_xdata(z)
    line2.set_ydata(fe)

    line21.set_xdata(z)
    line21.set_ydata(fe_osm)

    line22.set_xdata(z)
    line22.set_ydata(fe_surf)

    line3.set_xdata(z)
    line3.set_ydata(mobility_factor)

    line31.set_xdata(z)
    line31.set_ydata(conc_profile)

    ax[0].set_xlim(0, z_max)
    ax[0].set_ylim(0, phi_max)
    ax[1].set_xlim(0, z_max)
    ax[1].set_ylim(np.min([fe, fe_osm, fe_surf]), np.max([fe, fe_osm, fe_surf]))

    ax[2].set_ylim(0,np.max(conc_profile))


    fig.canvas.draw_idle()

chi_slider.on_changed(update)
chi_PC_slider.on_changed(update)
d_slider.on_changed(update)
ax[0].set_xlabel("$z$")
ax[0].set_ylabel(r"$\phi$")

ax[1].set_xlabel("$z$")
ax[1].set_ylabel(r"$F/k_B T$")

ax[2].set_xlabel("$z$")
ax[2].set_ylabel(r"$<c_{eq}>/c_{bulk}$")

ax[0].legend(loc = "lower left")
ax[1].legend(loc = "lower right")
# %%
