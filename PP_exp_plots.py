"""
@author: J.W.Spaak
Create the plots used in the manuscript
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np
import pandas as pd

import matplotlib as mpl
label_size = 18
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
plt.style.use('dark_background')

from NFD_for_experiments import NFD_experiment

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

def save_exp(ax_leg, handles, counter = [1]):
    print(counter[0])
    ax_leg.legend(handles = handles, fontsize = 16,
              bbox_to_anchor=(0.5,-0.1), loc="upper center", ncol = 5)
    plt.savefig("PP_slides/Experimental_data_{}.png".format(counter[0]),
                transparent = "True")
    counter[0] += 1
    
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

dens = np.array([BS4_dens[BS4_id, 2:inv_id], BS5_dens[~BS4_id, 2:inv_id]])

pars, N_t, N_t_data, fig, ax = NFD_experiment(dens, time[2:inv_id], r_i, 
                                              f0 = "linear")
N_star = pars["N_star"][[1,0],[0,1]]

# plot results
time_p = time_growth - time_growth[0] # start at day 0

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
###############################################################################
# star the actual plotting

# plot figure
fig, ax = plt.subplots(1,2,sharex = True, sharey = False, figsize = (11,9))

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

# axis ticks
ax[0].set_yticks([10**6, 10**7, 10**8, 10**9, N_star[0]])
ax[0].set_yticklabels([r"$10^6$", r"$10^7$", r"$10^8$", r"$10^9$", r"$N_1^*$"])
ax[0].set_xticks([0,20,40,60])
ax[1].set_xticks([0,20,40,60])
# increase size of N*_1
tick = ax[0].yaxis.get_major_ticks()[-1]
tick.label.set_fontsize(fs)
ax[1].set_yticks([N_star[1]])
ax[1].set_yticklabels([r"$N_2^*$"], fontsize = fs)

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

ax[1].set_axis_off()
ax[1].set_frame_on(False)
ax_arr1.set_axis_off()
ax_arr1.set_frame_on(False)
ax_leg = fig.add_subplot(111, frameon = False)
ax_leg.axis("off")

###############################################################################
# prepare axis
ax[1].set_axis_off()

# exp1 sp1
ax[0].errorbar(time_p, growth[0], growth_err[0], fmt = "^",  color = "red",
                      ecolor = "lightgrey")
time_f = np.linspace(time_growth[0],time[inv_id], 100)
ax[0].plot(time_f-time_growth[0], N_t["spec0_low"](time_f),
           color = "red")
arrow([time_p[0],growth[0,0]], [time_p[1],growth[0,1]], ax_arr0, 0.5, [3,0],
    arrowprops = dict(facecolor='orange'))
exp1, = ax_leg.plot(0, np.nan, "^",  color = "white",
                      label = "Exp1", markersize = 14)
handles = [exp1]
# plot equilibrium
ax[0].axhline(N_star[0], color = "white", linestyle = "-")
save_exp(ax_leg, handles)

# exp1 sp2, and invasion of sp1
ax[1].errorbar(time_p, growth[1], growth_err[1],fmt = "^",  color = "blue",
  ecolor = "lightgrey")
ax[1].plot(time_f-time_growth[0], N_t["spec1_low"](time_f),
           color = "blue")
ax[1].errorbar(time[inv_id:] - time_growth[0], BS4_inv[inv_id:],
  BS4_inv_err[inv_id:], fmt = "s",
  color = "red", ecolor = "lightgrey")
arrow([time_p[inv_id-2],BS4_inv[inv_id]], [time_p[inv_id-1],BS4_inv[inv_id+1]],
      ax_arr1, 0.5, [2,-0.1], arrowprops = dict(facecolor='purple'))
exp3, = ax_leg.plot(0, np.nan, "s",  color = "white",
                      label = "Invasion", markersize = 14)
handles = [exp1, exp3]
ax[1].axvline(time[inv_id] - time_growth[0], color = "grey", linestyle = ":")
ax[1].axhline(N_star[1], color = "white", linestyle = "-")
ax[1].set_axis_on()
ax[1].set_frame_on(True)

save_exp(ax_leg, handles)

# get growth for no-niche
t = 37
arrow([t - time_growth[0],N_t["spec0_low"](t)],
    [t + 7 - time_growth[0],N_t["spec0_low"](t + 7)],
    ax_arr0, 0.4, [0,0.2],
    arrowprops = dict(facecolor='lightblue', ls = 'dashed'))
save_exp(ax_leg, handles)

# growth arrow for second species, invasion
ax[0].errorbar(time[inv_id:] - time_growth[0], BS5_inv[inv_id:],
  BS4_inv_err[inv_id:], fmt = "s",
  color = "blue", ecolor = "lightgrey")
arrow([time_p[inv_id-2],BS5_inv[inv_id]], [time_p[inv_id-1],BS5_inv[inv_id+1]],
      ax_arr0, 0.2, [2,-0.1], arrowprops = dict(facecolor='purple'))
ax[0].axvline(time[inv_id] - time_growth[0], color = "grey", linestyle = ":")
arrow([time_p[0],growth[1,0]], [time_p[1],growth[1,1]], ax_arr1, 0.4, [2,0],
    arrowprops = dict(facecolor='orange'))
save_exp(ax_leg, handles)


# exp2 spec 1 and spec2
ax[0].errorbar(time_p, decline[0], decline_err[0], fmt = "o", color = "red"
  , ecolor = "lightgrey")
ax[0].plot(time_f-time_growth[0], N_t["spec0_high"](time_f), '--',
           color = "red")
ax[1].errorbar(time_p, decline[1], decline_err[1], fmt = "o", color = "blue"
  , ecolor = "lightgrey")
ax[1].plot(time_f-time_growth[0], N_t["spec1_high"](time_f), '--',
           color = "blue")
exp2, = ax_leg.plot(0, np.nan, "o",  color = "white",
                      label = "Exp2", markersize = 14)
handles = [exp1,exp3,exp2]
t = 44
arrow([t - time_growth[0],N_t["spec1_high"](t)],
    [t + 7 - time_growth[0],N_t["spec1_high"](t + 7)],
    ax_arr1, 0, [0,-0.2],
    arrowprops = dict(facecolor='lightblue', ls = 'dashed'))
save_exp(ax_leg, handles)
