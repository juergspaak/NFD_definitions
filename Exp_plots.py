"""
@author: J.W.Spaak
Create the plots used in the manuscript
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams["font.family"] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams.update({'font.size': 14})
import numpy as np
import pandas as pd

from NFD_for_experiments import NFD_experiment

file = "Exp_NFD,densities{}.csv"
BS4_dens = pd.read_csv(file.format("BS4"))
BS5_dens = pd.read_csv(file.format("BS5"))

BS4_dens = 1e3*BS4_dens.iloc[:,1:].values # convert to densities per ml
BS5_dens = 1e3*BS5_dens.iloc[:,1:].values # remove name column

# remove bad measurements
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

cut = 2 # cutoff first two datapoints, as species decrease

BS4_growth = geo_mean(BS4_dens[:3,cut:], axis = 0)
BS5_growth = geo_mean(BS5_dens[3:6,cut:], axis = 0)

BS4_decline = geo_mean(BS4_dens[6:9, cut:], axis = 0)
BS5_decline = geo_mean(BS5_dens[9:, cut:], axis = 0)

growth = np.array([BS4_growth, BS5_growth])
growth_err = np.array([np.nanstd(BS4_dens[:3,cut:], axis = 0),
                       np.nanstd(BS5_dens[3:6,cut:], axis = 0)])
decline = np.array([BS4_decline, BS5_decline])
decline_err = np.array([np.nanstd(BS4_dens[6:9,cut:], axis = 0),
                       np.nanstd(BS5_dens[9:12,cut:], axis = 0)])

# cut away data from equilibria
time_growth = time[cut:]
time_decline = time[cut:]

# compute invasion growth rate
BS4_id = np.array(2*(3*[True] + 3*[False]))
BS5_id = ~BS4_id
BS4_inv = geo_mean(BS4_dens[BS5_id], axis = 0)
BS5_inv = geo_mean(BS5_dens[BS4_id], axis = 0)
BS4_inv_err = np.nanstd(BS5_dens[BS4_id], axis = 0)
BS5_inv_err = np.nanstd(BS4_dens[BS5_id], axis = 0)
r_i = np.array([np.log(BS4_dens[BS5_id,19]/BS4_dens[BS5_id,17]),
                np.log(BS5_dens[BS4_id,19]/BS5_dens[BS4_id,17])])/7
r_i = np.nanmean(r_i, axis = -1)

dens = np.array([BS4_dens[BS4_id, cut:inv_id], BS5_dens[BS5_id, cut:inv_id]])
dil = np.log(0.9)/3.5 # dilution rate
# prepare
dens_flat = [4e10, 5e10]
gr = np.array([[dens_flat,len(dens_flat)*[dil]],
               [dens_flat,len(dens_flat)*[dil]]])

pars, N_t, N_t_data, fig, ax = NFD_experiment(dens, time[cut:inv_id], r_i, 
                                              f0 = "linear", s = "fac=0.95",
                                              growth_data = gr)
# compute R square
N_t_data = np.log(np.array([[N_t_data["spec0_low"], N_t_data["spec0_high"]],
                     [N_t_data["spec1_low"], N_t_data["spec1_high"]]]))
N_t_data.shape = 2,2,1,-1
# compute R2 in log:
dens2 = np.log(np.reshape(dens, (2,2,3,-1)))
axis = (2,3)
R2 = np.nansum((N_t_data-dens2)**2, axis = (axis))
R2 = 1-R2/np.nansum((dens2-np.nanmean(dens2, axis = (axis),
                                      keepdims = True))**2)

N_star = pars["N_star"][[1,0],[0,1]]


fig.savefig("NFD_computation_experiment.pdf")

# plot results
time_p = time_growth - time_growth[0] # start at day 0

###############################################################################
# start actual figure
# plot figure
fig, ax_all = plt.subplots(2,2,figsize = (9,9))
ax = ax_all[0]
b, w = 0.4, 0.45
#ax = [fig.add_axes([0, b, w, 1-b]),fig.add_axes([1-w, b, w, 1-b])]

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
ax[0].plot(time_f-time_growth[0], N_t["spec0_low"](time_f),
           color = "black", )
ax[0].plot(time_f-time_growth[0], N_t["spec0_high"](time_f), '--',
           color = "black")

ax[1].plot(time_f-time_growth[0], N_t["spec1_low"](time_f),
           color = "black")
ax[1].plot(time_f-time_growth[0], N_t["spec1_high"](time_f), '--',
           color = "black")


# plot invasion
ax[0].errorbar(time[inv_id:] - time_growth[0], BS5_inv[inv_id:],
  BS4_inv_err[inv_id:], fmt = "s",
  color = "grey", ecolor = "lightgrey")
ax[1].errorbar(time[inv_id:] - time_growth[0], BS4_inv[inv_id:],
  BS4_inv_err[inv_id:], fmt = "s",
  color = "black", ecolor = "lightgrey")

# customize axis labels etc
fs = 14 # fontsize
ax[0].set_ylim([8e5,6e9])
ax[0].semilogy()
ax[1].set_ylim([8e5,6e9])
ax[1].semilogy()
ax[0].set_xlim([min(time)-1,max(time)+1 -time_growth[0]+8])
ax[1].set_xlim([min(time)-1,max(time)+1 -time_growth[0]+8])
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
              bbox_to_anchor=(0.5,-0.08), loc="upper center", ncol = 5,
              fontsize = 11)

# plot equilibrium
ax[0].axhline(N_star[0], color = "black", linestyle = "-")
ax[1].axhline(N_star[1], color = "black", linestyle = "-")

ax[0].set_yticks([10**6, 10**7, 10**8, 10**9, N_star[0]])
ax[0].set_yticklabels([r"$10^6$", r"$10^7$", r"$10^8$", r"$10^9$", r"$N_1^*$"])
ax[0].set_xticks([0,20,40,50,60])
ax[0].set_xticklabels([0,20,40,"Invasion",60], fontsize = 12)
ax[1].set_xticks([0,20,40,50,60])
ax[1].set_xticklabels([0,20,40,"Invasion",60], fontsize = 12)
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
    
def arrow2(xy_tail, xy_head, ax, length = 1.5, shift_xy = 0,
          arrowprops = None, center = False):
    xy_tail = np.array(xy_tail)
    xy_head = np.array(xy_head)
    xy_tail[1] = np.log(xy_tail[1])
    xy_head[1] = np.log(xy_head[1])
    distance = xy_head-xy_tail
    axis_scale = [ax.get_xlim()[1]-ax.get_xlim()[0],
                  ax.get_ylim()[1]-ax.get_ylim()[0]]
    xy_head = xy_tail + length*distance/np.sqrt(np.sum((distance/axis_scale)**2))
    if center:
        diff = xy_head-xy_tail
        xy_head -= diff/2
        xy_tail -= diff/2
    ax.annotate(" ", xy_head + shift_xy, xy_tail + shift_xy,
                arrowprops = arrowprops)

# intrinsic growth rate
arrow([time_p[0],growth[0,0]], [time_p[1],growth[0,1]], ax_arr0, 1.5, [3,0],
    arrowprops = dict(facecolor='white'))

col = "None"
ax_arr0.text(5, 15, r"$f_1(0,0)$",  fontsize = fs-1,
            bbox = dict(facecolor=col, alpha=0.5, edgecolor = "None")) 
arrow([time_p[0],growth[1,0]], [time_p[1],growth[1,1]], ax_arr1, 1.1, [2,0],
    arrowprops = dict(facecolor='white'))
ax_arr1.text(5, 14, r"$f_2(0,0)$",  fontsize = fs-1,
            bbox = dict(facecolor=col, alpha=0.5, edgecolor = "None"))

# invasion growth
arrow([time_p[inv_id-2],BS5_inv[inv_id]], [time_p[inv_id-1],BS5_inv[inv_id+1]],
      ax_arr0, 1, [2,-0.1], arrowprops = dict(facecolor='black'))
ax_arr0.text(50, 14.9, r"$f_2(0,N_1^*)$",  fontsize = fs-2,
            bbox = dict(facecolor=col, alpha=0.5, edgecolor = "None"))
arrow([time_p[inv_id-2],BS4_inv[inv_id]], [time_p[inv_id-1],BS4_inv[inv_id+1]],
      ax_arr1, 2, [2,-0.1], arrowprops = dict(facecolor='black'))
ax_arr1.text(50, 15.1, r"$f_1(0,N_2^*)$",  fontsize = fs-2,
            bbox = dict(facecolor=col, alpha=0.5, edgecolor = "None"))

# no_niche growth rate
N_c = (pars["N_star"]*pars["c"])[[0,1],[1,0]]
t = 48
arrow([t - time_growth[0],N_t["spec0_low"](t)],
    [t + 7 - time_growth[0],N_t["spec0_low"](t + 7)],
    ax_arr0, 0.4, [0,0.2],
    arrowprops = dict(facecolor='white', ls = 'dashed'))
ax_arr0.text(t-15,18, r"$f_1(c_2N_2^*,0)$",  fontsize = fs-2,
            bbox = dict(facecolor=col, alpha=0.5, edgecolor = "None"))
t = 48
arrow([t - time_growth[0],N_t["spec1_high"](t)],
    [t + 7 - time_growth[0],N_t["spec1_high"](t + 7)],
    ax_arr1, 0.4, [0,-0.2],
    arrowprops = dict(facecolor='white', ls = 'dashed'))
ax_arr1.text(t-10,20.5, r"$f_2(c_1N_1^*,0)$",  fontsize = fs-2,
            bbox = dict(facecolor=col, alpha=0.5, edgecolor = "None"))

###############################################################################
# add axes for resource use and for growth rates
# load absorption spectra
k_spec = pd.read_csv("Lights_and_absorptions.csv")
n_bars = 15
#ax_res = fig.add_axes([0,0,w, b-0.1])
ax_res = ax_all[1,0]
ax_res.set_title("C: Absorption spectrum", loc = "left")
ax_res.set_xlabel(r"wavelength $[nm]$")
ax_res.set_ylabel(r"absorption $[ml\cdot cells^{-1}\cdot m^{-1}]$")

dist = len(k_spec)//n_bars
# change units from absorption spectrum, current unit: mul/cells/cm
# change to ml/cells/m
unit_conv = 1000*100
ax_res.bar(k_spec["lambda"][::dist], unit_conv*(k_spec["BS4; 0"])[::dist],
           width = dist,
           color = "black")
ax_res.bar(k_spec["lambda"][::dist]+dist, unit_conv*(k_spec["BS5; 0"])[::dist],
           width = dist,
           color = "lightgrey")
ax_res.set_yticks([0,0.3,0.6])

ax_gro = ax_all[1,1]
locations = np.arange(5)*3
# corresponds to c_0
c_estimated = np.sqrt((np.sum(k_spec["BS4; 0"]**2*k_spec["I_in"])
                /np.sum(k_spec["BS5; 0"]**2*k_spec["I_in"])))


ax_gro.barh(locations+1, [pars["f0"][0], pars["r_i"][0],
                       pars["fc"][0], pars["c"][1,0], 1/c_estimated],
     color = "black")

ax_gro.barh(locations, [pars["f0"][1], pars["r_i"][1],
                       pars["fc"][1], pars["c"][0,1], c_estimated],
     color = "lightgrey")

values = np.array([pars["f0"], pars["r_i"],
                       pars["fc"], pars["c"][[1,0],[0,1]],
                       [1/c_estimated, c_estimated]]).T
for i,rect in enumerate(ax_gro.patches):
    if not(i%5==2):
        continue
    ax_gro.text(0.3, rect.get_y()+rect.get_width()/2,
                np.round(values.reshape(-1)[i],3), ha = "right")

names = [r"$f_i(0,0)$", r"$f_i(0,N_j^*)$", r"$f_i(c_jN_j^*,0)$", r"$c$", 
         "tot.\nconsumption"]
names = ["intrinsic\ngrowth", "invasion\ngrowth", "no-niche\ngrowth",
         r"$c$", "relative total\nabsorption"]
lines = np.linspace(*ax_gro.get_ylim(), len(names)+1)
for line in lines:
    ax_gro.axhline(line, color = "k", linestyle = ":")

ax_gro.set_xticks([0, 0.5, 1.0])
#ax_gro.set_yticklabels(names, fontsize = 12)
name_positions = np.linspace(*ax_gro.get_ylim(), 2*len(names)+1)
ax_gro.set_yticks(name_positions[1::2])
ax_gro.set_yticklabels(names, fontsize = 12)
ax_gro.set_title("D: growth rates", loc = "left")
fig.tight_layout()
fig.savefig("Experimental_data.pdf")
fig.savefig("Experimental_data.eps")
#"""
plt.show()
print(pars["ND"])
print(pars["FD"])
