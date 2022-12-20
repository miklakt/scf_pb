#%%
import scf_pb
import numpy as np
import time
import numpy as np
import matplotlib.pyplot as plt
#matplotlib settings
LAST_USED_COLOR = lambda: plt.gca().lines[-1].get_color()
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
import matplotlib.pyplot as plt
import itertools
import matplotlib.style as style
style.use('tableau-colorblind10')
mpl_markers = ('o', '+', 'x', 's', 'D')
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
N=2000
sigma = 0.08

a0=0.18
a1=-0.09
d=np.arange(1, 31, step=1)

chi_major = [0.00, 0.25, 0.50, 0.75, 1.00]
chi_ticks= list(np.arange(0,1.05, 0.05))
#chi = np.sort(list(set(chi_major+chi_ticks)))
chi = chi_ticks
chi_crit = scf_pb.chi_PC_critical_v(
    N=N, sigma = sigma, 
    chi = chi,
    a0=a0, a1=a1, 
    d=d, 
    k_smooth = 1, 
    a=0, b = "H",
    progressbar = True,
    )
#%%
for chi_crit_ in chi_crit:
    cut_off_idx = np.argmax(chi_crit_)
    chi_crit_[:-(len(chi_crit_)-cut_off_idx)] = np.nan

#%%
def find_nearest_value_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_nearest_chi_PS(d, chi_PC, d_arr, chi_PS_arr, chi_crit):
    d_index = np.where(d_arr==d)
    chi_PS_on_chi_PC =chi_crit.T[d_index]
    chi_PS_idx = find_nearest_value_idx(chi_PS_on_chi_PC, chi_PC)
    return chi_PS_arr[chi_PS_idx]

# %%
fig, ax  = plt.subplots()
markers = itertools.cycle(mpl_markers)
for chi_crit_, chi_ps in zip(chi_crit, chi):
    #d_ = d[cut_off_idx:]
    if chi_ps in chi_major:
            ax.plot(
            d, 
            chi_crit_,
            label = f"{chi_ps:.2f}", 
            linewidth= 2,
            marker = next(markers),
            markerfacecolor = "None",
            markersize = 8,
            markeredgewidth=2,
            markevery = 0.2,
            )
    else:
        ax.plot(d, chi_crit_,
        color = "black", linewidth = 0.2)

#[ax.plot(chi_crit_,d, color = "black") for chi_crit_, d in zip(chi_crit.T, d)] 

ax.legend(title = r"$\chi_{\rm PS}$")
ax.grid()
ax.set_ylabel(r"$\chi_{\rm PC}$")
ax.set_xlabel(r"$d$")
ax.set_ylim(-3, -1.2)

fig.set_size_inches(3,4)

fig.savefig("fig8a.svg", bbox_inches = "tight")
fig.savefig("fig8a.pdf", bbox_inches = "tight")


#%%
fig, ax  = plt.subplots()
#markers = itertools.cycle(mpl_markers)
for chi_crit_, d_ in zip(chi_crit.T, d):
    if d_ in [2.0, 4.0, 8.0,  16.0,  30.0]:
        ax.plot( chi, chi_crit_, label = d_, linewidth = 2, #linestyle = "--",
            marker = next(markers),
            #markerfacecolor = "None",
            markersize = 8,
            #markeredgewidth=2,
            #markevery = 0.2,
            )
    else:
        ax.plot( chi, chi_crit_, color = 'black', linewidth = 0.2)

#[ax.plot(chi_crit_,d, color = "black") for chi_crit_, d in zip(chi_crit.T, d)] 

ax.legend(title = r"d", loc = 'lower left')
ax.grid()
#ax.set_ylabel(r"$\chi_{\rm PC}$")
ax.set_xlabel(r"$\chi_{\rm PS}$")
ax.set_ylim(-3, -1.2)
fig.set_size_inches(3,4)

fig.savefig("fig8b.svg", bbox_inches = "tight")
fig.savefig("fig8b.pdf", bbox_inches = "tight")
# %%
