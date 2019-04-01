"""
@author: J.W.Spaak
Create the plots used in the manuscript
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import rainbow

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


# fill gaps and convolve
n = 5 #number of running average
for dens in [BS4_growth, BS5_growth, BS4_decline, BS5_decline]:
    dens[2:-2] = np.convolve(dens, np.ones((n,))/n, mode='same')[2:-2]

growth = np.array([BS4_growth, BS5_growth])
decline = np.array([BS4_decline, BS5_decline])

# cut away data from equilibria
time_growth = time[2:]
time_decline = time[2:]

# compute invasion growth rate
BS4_id = np.array(2*(3*[True] + 3*[False]))
BS5_id = ~BS4_id
r_i = np.array([np.log(BS4_dens[BS5_id,19]/BS4_dens[BS5_id,17]),
                np.log(BS5_dens[BS4_id,19]/BS5_dens[BS4_id,17])])/6
r_i = np.nanmean(r_i, axis = -1)
data_growth = np.arange(0,14)
data_decline = np.arange(0,14)
N_star = (growth[:,data_growth[-1]] + decline[:,data_decline[-1]])/2

# plot figure
fig, ax = plt.subplots(1,2,sharex = True, sharey = True, figsize = (11,9))

ax[0].semilogy(time, BS4_dens[0:3].T, 'go')
ax[0].semilogy(time, BS4_dens[6:9].T, "g^")
ax[0].semilogy(time[inv_id:], BS5_dens[BS4_id,inv_id:].T, 'rs')
#ax[0].semilogy(time_growth, BS4_growth, 'o', color = "black")
#ax[0].semilogy(time_decline, BS4_decline, '^', color = "black")

ax[1].semilogy(time, BS5_dens[3:6].T, 'ro')
ax[1].semilogy(time, BS5_dens[9:12].T, "r^")
ax[1].semilogy(time[inv_id:], BS4_dens[BS5_id,inv_id:].T, 'gs')
#ax[1].semilogy(time_growth, BS5_growth, 'o', color = "black")
#ax[1].semilogy(time_decline, BS5_decline, '^', color = "black")

# add division into different experiments
fs = 16 # fontsize
ax[0].set_ylim([8e5,6e9])
ax[0].set_xlim([min(time)-1,max(time)+1])
ax[0].hlines(N_star[0],*ax[0].get_xlim(),
      color = "black", linestyles = "dotted")
ax[0].vlines(time[inv_id], *ax[0].get_ylim(), color = "black",
      linestyles = "dotted")

ax[0].text(10,1.5e6,"   Monoculture growth\n   from low abundance", ha = "left", 
  va = "center", fontsize = fs, bbox = {"fc":"None", "ec":"black"})
ax[0].plot(12,1.5e6, 'go', markersize = 10)
ax[0].text(5,4e9,"   Monoculture growth\n   from high abundance", va = "center",
  ha = "left", fontsize = fs, bbox = {"fc":"None", "ec":"black"})
ax[0].plot(7,4e9, 'g^', markersize = 10)
ax[0].text(time[inv_id]-4,np.nanmean(BS5_dens[BS4_id,inv_id]),
  "    Invasion", va = "center", ha = "center", fontsize = fs,rotation = 90,
  bbox = {"fc":"None", "ec":"black"})
ax[0].plot(time[inv_id]-4.5,np.nanmean(BS5_dens[BS4_id,inv_id])/1.8,
  'rs', markersize = 10)
plt.yticks([1e6,1e7,1e8,N_star[0], 1e9],
          [r"$10^6$", r"$10^7$", r"$10^8$", r"$N^*$",r"$10^9$"], fontsize = fs)


ax[1].hlines(N_star[1],*ax[1].get_xlim(),
      color = "black", linestyles = "dotted")
ax[1].vlines(time[inv_id], *ax[1].get_ylim(), color = "black",
      linestyles = "dotted")

ax[0].set_xlabel("Time [days]", fontsize = fs)
ax[1].set_xlabel("Time [days]", fontsize = fs)
ax[0].set_ylabel("Density [cells/ml]", fontsize = fs)

ax[0].set_title("A")
ax[1].set_title("B")

fig.savefig("Experimental_data.pdf")


# compute NFD values
pars, fig, ax = NFD_experiment(N_star, time_growth[data_growth], 
            growth[:,data_growth], time_decline[data_decline],
            decline[:,data_decline], r_i, k = 3)

ax[0,0].semilogy()
ax[0,1].semilogy()

ax[1,0].semilogx()
ax[1,1].semilogx()

colors = rainbow(np.linspace(0,1,len(BS4_growth)))

for i in range(len(BS4_growth)):
    ax[1,0].scatter(BS4_growth,[pars["f"]([N,0])[0] for N in BS4_growth],
                    c = np.linspace(0,1,len(BS4_growth)))
    ax[1,0].scatter(BS4_decline,[pars["f"]([N,0])[0] for N in BS4_decline],
                    c = np.linspace(0,1,len(BS4_decline)), marker = '^')
    
    ax[1,1].scatter(BS5_growth,[pars["f"]([0,N])[1] for N in BS5_growth],
                    c = np.linspace(0,1,len(BS4_growth)))
    ax[1,1].scatter(BS5_decline,[pars["f"]([0,N])[1] for N in BS5_decline],
                    c = np.linspace(0,1,len(BS5_decline)), marker = '^')
fig.savefig("NFD_computation_experiment.png")

