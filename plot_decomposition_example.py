import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations

import Letten_2018_functions as lf
from scheme_plot import plot_scheme

n_spec = 4
combs = []
# create all species combinations
for i in range(n_spec):
    combs.extend(combinations(np.arange(n_spec), i))
    
combs = [list(x) for x in combs]
combs.append(np.arange(4))


r_i = np.empty((len(combs), n_spec))
dens_all = []

for i, comb in enumerate(combs):
    print(i)
    dens_precise, time_precise = lf.get_equilibrium_densities(comb)
    dens_all.append(dens_precise)
    

    r_i[i] = lf.decomposition(dens_precise, time_precise)["real"]

dens_all = np.array(dens_all)

#"""

fig, ax = plt.subplots(4,4, sharex = True, figsize = (10,10))
ax = ax.flatten()
for i in range(len(combs)):
    ax[i].plot(time_precise, dens_all[i].T, label = "hi")
    ax[i].set_title(str(combs[i]))
ax[0].legend(labels = ["Sp {}".format(i) for i in range(n_spec)] + ["Res"])

fig.tight_layout()#"""

# compute invasion scheme
scheme = r_i.copy()
mean_dens = np.mean(dens_all, axis = -1)
for i in range(len(combs)):
    if (np.mean(dens_all[i, combs[i]], axis = -1)<1).any():
        scheme[i] = np.nan # remove subcommunity
    else:
        scheme[i, combs[i]] = 0 # set invasion growth to actual zero
        
scheme = scheme[np.isfinite(scheme[:,0])]



# compute growth rates decompositions
decomps = np.empty((3,2,6))
keys = ["const", "case", "res", "rand", "real"]

def decomp(r_is, inv, res):
    epsilon = [r_is["const"],
                    r_is["case"] - r_is["const"],
                    r_is["res"] - r_is["const"],
                    r_is["rand"] - r_is["case"] - r_is["res"] + r_is["const"],
                    0,0]
    epsilon[-2] = r_is["real"] - r_is["rand"]
    epsilon[-1] = r_is["real"]
    epsilon = np.array(epsilon)
    return (epsilon[:, inv] - np.mean(epsilon[:, res], axis = 1, keepdims = True)).T

# decomposition for species 2
inv = [2,3]
res = [0]
i = 1 # corresponds to community where only species 1 is present
r_is = lf.decomposition(dens_all[i], time_precise)
decomps[0] = decomp(r_is, inv, res)

# decomposition for species 1
inv = [2,3]
res = [1]
i = 2 # corresponds to community where only species 2 is present
r_is = lf.decomposition(dens_all[i], time_precise)
decomps[2] = decomp(r_is, inv, res)

# decomposition for species 1
inv = [2,3]
res = [0,1]
i = 5 # corresponds to community where only species 2 is present
r_is = lf.decomposition(dens_all[i], time_precise)
decomps[1] = decomp(r_is, inv, res)


##############################################################################•
# plot results
fig, ax = plt.subplots(2,2, figsize = (9,9), sharex = "row")


x = np.arange(decomps.shape[-1])[::-1]
h = 0.3
add = [-0.2,0,0.2]
label = ["Sp. 1 resident",
         "Sp. 1 & 2 resident",
         "Sp. 2 resident"]
for i in range(3):
    ax[1,0].barh(x + add[i], decomps[i,0], height = h, label = label[i])
    ax[1,1].barh(x + add[i], decomps[i,1], height = h)


"""ax[1,1].barh(x, decomps[0], height = h)
ax[1,1].barh(x+0.2, decomps[1], height = h)

ax[1,0].barh(x, decomps[2], height = h)
ax[1,0].barh(x+0.2, decomps[3], height = h)"""

ax[1,0].set_yticks(x)
ax[1,0].set_yticklabels([r"$r_i$",
                         r"$\Delta_i^{(RS)}$",
                         "$\Delta_i^{(R{\#}S)}$",
                         r"$\Delta_i^{R}$",
                         r"$\Delta_i^{S}$",
                         r"$\Delta_i^{0}$"][::-1], fontsize = 16)
ax[1,1].set_yticks(x)
ax[1,1].set_yticklabels([])
ax[1,0].set_xlim([-0.07,0.05])
ax[1,0].legend(loc = "center left")
ax[1,0].set_xlabel("Growth rate")
ax[1,1].set_xlabel("Growth rate")
ax[0,0].set_xticks([])
ax[0,0].set_ylabel("Species richness")
ax[1,0].set_title("Species 3")
ax[1,1].set_title("Species 4")
ax[0,0].set_title("Full invasion graph")
ax[0,1].set_title("Invasion graph\nof permanent sub-community")

plot_scheme(scheme, ax[0,0], sp_names = "index_at_1")
plot_scheme(scheme[[0,1,2,-2], :2],ax[0,1], sp_names = "index_at_1")

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCD"[i], loc = "left")

fig.savefig("Figure_decomposition_example.pdf")