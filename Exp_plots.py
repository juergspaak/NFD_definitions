"""
@author: J.W.Spaak
Create the plots used in the manuscript
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams["font.family"] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
import numpy as np
import pandas as pd

from NFD_for_experiments import NFD_experiment

file = "Exp_NFD,densities{}.csv"
BS4_dens = pd.read_csv(file.format("BS4"))
BS5_dens = pd.read_csv(file.format("BS5"))

BS4_dens = 1e3*BS4_dens.iloc[:,1:].values # convert to densities per ml
BS5_dens = 1e3*BS5_dens.iloc[:,1:].values

BS4_dens[7, 9:15] = np.nan
BS4_dens[8,15:] = np.nan

BS4_dens[1, 17] = np.nan
BS4_dens[5, 18] = np.nan
BS5_dens[1, 17] = np.nan
BS4_dens[2, 19] = np.nan

# measurment taken all 3.5 days
time = 3.5*np.arange(BS4_dens.shape[1]-1)
# 1. measurment after the invasion taken 1 day after the invasion
inv_id = 17 # invasion occured in the 18 measurment
time = np.insert(time, inv_id, (inv_id - 1)*3.5 + 1)

def geo_mean(data, axis = None):
    return np.nanmean(data, axis = axis)
    return np.exp(np.nanmean(np.log(data), axis = axis))


BS4_growth = geo_mean(BS4_dens[:3,2:], axis = 0)
BS5_growth = geo_mean(BS5_dens[3:6,2:], axis = 0)

BS4_decline = geo_mean(BS4_dens[6:9, 2:], axis = 0)
BS5_decline = geo_mean(BS5_dens[9:, 2:], axis = 0)

growth = np.array([BS4_growth, BS5_growth])
growth_err = np.array([np.nanstd(BS4_dens[:3,2:], axis = 0),
                       np.nanstd(BS5_dens[3:6,2:], axis = 0)])
decline = np.array([BS4_decline, BS5_decline])
decline_err = np.array([np.nanstd(BS4_dens[6:9,2:], axis = 0),
                       np.nanstd(BS5_dens[9:12,2:], axis = 0)])

# cut away data from equilibria
time_growth = time[2:]
time_decline = time[2:]

# compute invasion growth rate
BS4_id = np.array(2*(3*[True] + 3*[False]))
BS5_id = ~BS4_id
BS4_inv = geo_mean(BS4_dens[BS5_id], axis = 0)
BS5_inv = geo_mean(BS5_dens[BS4_id], axis = 0)
BS4_inv_err = np.nanstd(BS5_dens[BS4_id], axis = 0)
BS5_inv_err = np.nanstd(BS4_dens[BS5_id], axis = 0)
r_i = np.array([np.log(BS4_dens[BS5_id,19]/BS4_dens[BS5_id,17]),
                np.log(BS5_dens[BS4_id,19]/BS5_dens[BS4_id,17])])/6
r_i = np.nanmean(r_i, axis = -1)
data_growth = np.arange(0,16)
data_decline = np.arange(0,16)
N_star = (growth[:,data_growth[-1]] + decline[:,data_decline[-1]])/2

# compute NFD values
pars, N_t, fig, ax = NFD_experiment(N_star, time_growth[data_growth], 
            growth[:,data_growth], time_decline[data_decline],
            decline[:,data_decline], r_i)

fig.savefig("NFD_computation_experiment.pdf")

# plot results
time_p = time_growth - time_growth[0] # start at day 0


# plot figure
fig, ax = plt.subplots(1,2,sharex = True, sharey = False, figsize = (11,9))

# plot the experimentaly measured data
ax[0].errorbar(time_p, growth[0], growth_err[0], fmt = "^",  color = "black",
                      ecolor = "lightgrey")

ax[0].errorbar(time_p, decline[0], decline_err[0], fmt = "o", color = "black"
  , ecolor = "lightgrey")

ax[1].errorbar(time_p, growth[1], growth_err[1],fmt = "^",  color = "grey",
  ecolor = "lightgrey")
ax[1].errorbar(time_p, decline[1], decline_err[1], fmt = "o", color = "grey"
  , ecolor = "lightgrey")

# plot the fitted data
time_f = np.linspace(time_growth[0],time[inv_id], 100)
ax[0].plot(time_f-time_growth[0], [N_t["exp1_spec0"](t) for t in time_f],
           color = "black", )
ax[0].plot(time_f-time_growth[0], [N_t["exp2_spec0"](t) for t in time_f], '--',
           color = "black")

ax[1].plot(time_f-time_growth[0], [N_t["exp1_spec1"](t) for t in time_f],
           color = "black")
ax[1].plot(time_f-time_growth[0], [N_t["exp2_spec1"](t) for t in time_f], '--',
           color = "black")


# plot invasion
ax[0].errorbar(time[inv_id:] - time_growth[0], BS5_inv[inv_id:],
  BS4_inv_err[inv_id:], fmt = "s",
  color = "grey", ecolor = "lightgrey")
ax[1].errorbar(time[inv_id:] - time_growth[0], BS4_inv[inv_id:],
  BS4_inv_err[inv_id:], fmt = "s",
  color = "black", ecolor = "lightgrey")

# customize axis labels etc
fs = 16 # fontsize
ax[0].set_ylim([8e5,6e9])
ax[0].semilogy()
ax[1].set_ylim([8e5,6e9])
ax[1].semilogy()
ax[0].set_xlim([min(time)-1,max(time)+1 -time_growth[0]+8])
ax[0].set_xlabel("Time [days]", fontsize = fs)
ax[1].set_xlabel("Time [days]", fontsize = fs)
ax[0].set_ylabel("Density [cells/ml]", fontsize = fs)

# add title and legend
ax[0].set_title("A; Species 1")
ax[1].set_title("B; Species 2")
ax_leg = fig.add_subplot(111, frameon = False)
ax_leg.axis("off")
exp1, = ax_leg.plot(0, np.nan, "-^",  color = "black", fillstyle = "none",
                      label = "Monoculture growth\nfrom low abundance")
exp2, = ax_leg.plot(0, np.nan, "--o",  color = "black", fillstyle = "none",
                      label = "Monoculture growth\nfrom high abundance")
exp3, = ax_leg.plot(0, np.nan, "s",  color = "black", fillstyle = "none",
                      label = "Invasion")
spec1 = mpatches.Patch(color='black', label = 'Species 1')
spec2 = mpatches.Patch(color='grey', label = 'Species 2')
ax_leg.legend(handles = [exp1,exp2,exp3,spec1, spec2],
              bbox_to_anchor=(0.5,-0.1), loc="upper center", ncol = 5)

# plot equilibrium
ax[0].axhline(N_star[0], color = "black", linestyle = "-")
ax[1].axhline(N_star[1], color = "black", linestyle = "-")

ax[0].set_yticks([10**6, 10**7, 10**8, 10**9, N_star[0]])
ax[0].set_yticklabels([r"$10^6$", r"$10^7$", r"$10^8$", r"$10^9$", r"$N_1^*$"])
ax[0].set_xticks([0,20,40,50,60])
ax[0].set_xticklabels([0,20,40,"Invasion",60])
ax[1].set_xticks([0,20,40,50,60])
ax[1].set_xticklabels([0,20,40,"Invasion",60])
# increase size of N*_1
tick = ax[0].yaxis.get_major_ticks()[-1]
tick.label.set_fontsize(fs)
ax[1].set_yticks([N_star[1]])
ax[1].set_yticklabels([r"$N_2^*$"], fontsize = fs)

# plot invasion
ax[0].axvline(time[inv_id] - time_growth[0], color = "grey", linestyle = ":")
ax[1].axvline(time[inv_id] - time_growth[0], color = "grey", linestyle = ":")

# add mathematical expressions
ax_arr0 = ax[0].twinx()
ax_arr0.set_ylim(*np.log(ax[0].get_ylim()))
ax_arr0.set_xlim(ax[0].get_xlim())
ax_arr0.set_yticks([])
ax_arr1 = ax[1].twinx()
ax_arr1.set_ylim(*np.log(ax[1].get_ylim()))
ax_arr1.set_xlim(ax[1].get_xlim())
ax_arr1.set_yticks([])
gr_arr = np.log(growth)

# to draw arrows
def arrow(xy_tail, xy_head, ax, increase_length = 0, shift_xy = 0,
          arrowprops = None):
    xy_tail = np.array(xy_tail)
    xy_head = np.array(xy_head)
    xy_tail[1] = np.log(xy_tail[1])
    xy_head[1] = np.log(xy_head[1])
    xy_tail = xy_tail + (xy_tail-xy_head)*increase_length
    ax.annotate(" ", xy_head + shift_xy, xy_tail + shift_xy,
                arrowprops = arrowprops)

# no competition growth
arrow([time_p[0],growth[0,0]], [time_p[1],growth[0,1]], ax_arr0, 0.5, [3,0],
    arrowprops = dict(facecolor='white'))
ax_arr0.text(5, 15, r"$f_1(0,0)$",  fontsize = fs-1,
            bbox = dict(facecolor='white', alpha=0.5, edgecolor = "None"))    
arrow([time_p[0],growth[1,0]], [time_p[1],growth[1,1]], ax_arr1, 0.4, [2,0],
    arrowprops = dict(facecolor='white'))
ax_arr1.text(5, 14, r"$f_2(0,0)$",  fontsize = fs-1,
            bbox = dict(facecolor='white', alpha=0.5, edgecolor = "None"))

# invasion growth
arrow([time_p[inv_id-2],BS5_inv[inv_id]], [time_p[inv_id-1],BS5_inv[inv_id+1]],
      ax_arr0, 0.2, [2,-0.1], arrowprops = dict(facecolor='black'))
ax_arr0.text(50, 15.9, r"$f_2(0,N_1^*)$",  fontsize = fs-2,
            bbox = dict(facecolor='white', alpha=0.5, edgecolor = "None"))
arrow([time_p[inv_id-2],BS4_inv[inv_id]], [time_p[inv_id-1],BS4_inv[inv_id+1]],
      ax_arr1, 0.5, [2,-0.1], arrowprops = dict(facecolor='black'))
ax_arr1.text(50, 16.1, r"$f_1(0,N_2^*)$",  fontsize = fs-2,
            bbox = dict(facecolor='white', alpha=0.5, edgecolor = "None"))

# no_niche growth rate
t = 37
arrow([t - time_growth[0],N_t["exp1_spec0"](t)],
    [t + 7 - time_growth[0],N_t["exp1_spec0"](t + 7)],
    ax_arr0, 0.4, [0,0.2],
    arrowprops = dict(facecolor='white', ls = 'dashed'))
ax_arr0.text(t-15,19.4, r"$f_1(c_2N_2^*,0)$",  fontsize = fs-2,
            bbox = dict(facecolor='white', alpha=0.5, edgecolor = "None"))
t = 44
arrow([t - time_growth[0],N_t["exp2_spec1"](t)],
    [t + 7 - time_growth[0],N_t["exp2_spec1"](t + 7)],
    ax_arr1, 0, [0,-0.2],
    arrowprops = dict(facecolor='white', ls = 'dashed'))
ax_arr1.text(t-15,19.6, r"$f_2(c_1N_1^*,0)$",  fontsize = fs-2,
            bbox = dict(facecolor='white', alpha=0.5, edgecolor = "None"))

fig.savefig("Experimental_data.pdf", bbox_inches = "tight")
