"""
@author: J.W.Spaak
Create Fiugre 5
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams["font.family"] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams.update({'font.size': 14})
import numpy as np

from numerical_NFD import NFD_model
from scipy.integrate import odeint

def arrow(xy_tail, xy_head, ax, length = 1.5, shift_xy = 0,
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
c_help = 0.5
A = np.array([[0.25*c_help, 0.3], [0.3125*c_help, 1.5]])/2
mu = 1

def LV(N, t = 0):
    return 1-A.dot(N)

def LV_ode(N,t):
    return N*LV(N,t)

pars = NFD_model(LV)

fig, ax = plt.subplots(1,2,figsize = (9,5), sharex = "row")
ax[0].set_ylim([0.01, 40])
ax[1].set_ylim([0.01, 40])
ax[0].set_xlim([-0.5,11])
fs = 14
colors = ["black", "grey"]


for i in range(2):
    time = np.linspace(0,8,300)
    step_t = 20
    j = 1-i
    dens = np.zeros(2)
    # plot monoculture from low densities
    dens[i] = 0.1
    growth = odeint(LV_ode, dens, time)[:,i]
    ax[i].plot(time, growth, linestyle = "--",
      color = colors[i])
    # plot monoculture from high densities
    dens[i] = 30
    decline = odeint(LV_ode, dens, time)[:,i]
    ax[i].plot(time, decline,
      color = colors[i], linestyle = ':')
    ax[i].axhline(pars["N_star"][j,i], color = "black")
    
    # plot the invasion growth rate
    time_inv = np.linspace(time[-1],time[-1]*1.2,3)
    dens[i] = pars["N_star"][j,i]
    dens[j] = 0.1
    invasion = odeint(LV_ode, dens, time_inv)[:,j]
    ax[i].plot(time_inv, invasion,
      color = colors[j], linestyle = '-')
    ax[i].semilogy()
    
    ax[i].axvline(min(time_inv), linestyle = ':', color = "black")
    
    ax_arr = ax[i].twinx()
    ax_arr.set_ylim(*np.log(ax[i].get_ylim()))
    ax_arr.set_xlim(ax[i].get_xlim())
    ax_arr.set_yticks([])
    
    # intrinsic growth rates
    arrow([time[0],growth[0]], [time[1],growth[1]], ax_arr, 0.2, [0,0],
          arrowprops = dict(facecolor='white'))
    ax_arr.text(time[0]+0.8, np.log(growth[0])+0.2, r"$f_{}(0,0)$".format(i+1), 
            fontsize = fs-1,
            bbox = dict(facecolor='white', alpha=0.0, edgecolor = "None"))  
    
    # invasion growth rates
    arrow([time_inv[0],invasion[0]], [time_inv[1],invasion[1]], ax_arr, 0.2,
          [0,0.2],
          arrowprops = dict(facecolor='black'))
    ax_arr.text(time_inv[0], np.log(invasion[0]) - (-1)**(i+1)*0.3,
                r"$f_{}(0,N_{}^*)$".format(j+1,i+1), fontsize = fs-1,
            bbox = dict(facecolor='white', alpha=0.0, edgecolor = "None"))
    
    # no-niche growth rate
    cjNj = (pars["N_star"]*pars["c"])[i,j]
    if cjNj>pars["N_star"][j,i]:
        t_id = np.argmax(decline<cjNj)
        no_niche = decline
        step = 1
    else:
        t_id = np.argmax(growth>cjNj)
        no_niche = growth
        step = -1
    
    # no-niche growth rate
    arrow([time[t_id],no_niche[t_id]], [time[t_id+1],no_niche[t_id+1]], ax_arr,
          0.2, [0,0.0], center = True,
          arrowprops = dict(facecolor='white', ls = 'dashed'))
    ax_arr.text(time[t_id]+0.8, np.log(no_niche[t_id]) - step*0,
                r"$f_{}(c_{}N_{}^*,0)$".format(i+1,j+1,j+1), fontsize = fs-1,
            bbox = dict(facecolor='white', alpha=0.0, edgecolor = "None"))
    ax[i].set_yticks([pars["N_star"][j,i], cjNj])
    ax[i].set_yticklabels([r'$N_{}^*$'.format(i+1),
                              r'$c_{}N_{}^*$'.format(j+1,j+1)],
      fontsize = fs)

ax[0].set_ylabel("Density (log scale)")
ax[0].set_title("A: Species 1", loc = "left")
ax[1].set_title("B: Species 2", loc = "left")

ax[0].set_xticks([0,time[-1]])
ax[0].set_xticklabels([0,"Invasion"], fontsize = 12)
ax[1].set_xticks([0,time[-1]])
ax[1].set_xticklabels([0,"Invasion"], fontsize = 12)

###############################################################################
# legend
ax_leg = fig.add_subplot(111, frameon = False)
ax_leg.axis("off")
exp1, = ax_leg.plot(0, np.nan, "--",  color = "black", fillstyle = "none",
                      label = "Monoculture growth\nfrom low abundance")
exp2, = ax_leg.plot(0, np.nan, ":",  color = "black", fillstyle = "none",
                      label = "Monoculture growth\nfrom high abundance")
spec1 = mpatches.Patch(color='black', label = 'Species 1')
spec2 = mpatches.Patch(color='grey', label = 'Species 2')
ax_leg.legend(handles = [exp1,exp2,spec1, spec2]
              , loc = "center", bbox_to_anchor=(0.5,-0.17), ncol = 4,
              fontsize = 12)
fig.tight_layout()
fig.savefig("Figure5.pdf")