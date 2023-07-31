import numpy as np
import matplotlib.pyplot as plt



from nfd_definitions.numerical_NFD import NFD_model
import pandas as pd

import scheme_plot as schp
    

A = pd.read_csv("Geijzendorffer 2011.csv", index_col = 0).values
n = len(A)
sp_names = np.arange(len(A)) + 1


r_i, equi_all = schp.compute_scheme_for_LV(A)
ind = np.all(equi_all>=0, axis = 1)
equi_all = equi_all[ind]
pars = schp.compute_NFD_LV(A)

# compute niche and fitness difference at all equilibria
mu_i = np.ones(r_i.shape)
eta_i = mu_i - np.diag(A)*np.einsum("ij,nj->ni", pars["c"], equi_all)

ND_all = (r_i - eta_i)/(mu_i - eta_i)
F_all = -eta_i/(mu_i - eta_i)

ND_all[np.isclose(r_i, 0)] = np.nan
F_all[np.isclose(r_i, 0)] = np.nan

richness = np.sum(np.isclose(r_i, 0), axis = 1)

fig, ax = plt.subplots(2,3, sharex = True, sharey = True,
                       figsize = (9,9))

for a in ax[-1]:
    a.set_xlabel("Niche differences")
for a in ax[:,0]:
    a.set_ylabel("Fitness differnecs")
ax = ax.flatten()

titles = ["Richness {}".format(i) for i in range(1,6)] +["All Richness"]
for i,a in enumerate(ax):
    a.set_title("ABCDEF"[i], loc = "left")
    a.plot([0,1],[0,1], 'r')
    a.axhline(0, color = "grey")
    a.axvline(0, color = "grey")
    a.axhline(1, color = "grey")
    a.axvline(1, color = "grey")
    a.set_title(titles[i])

text_args = dict(ha = "center", va = "center", fontsize = 12)
cmap = plt.get_cmap("viridis")
colors = cmap(np.linspace(0,1, len(A)))

# remove NFD values outside of the figure boundaries
ND_lim = [-0.1, 2.5]
ND_all[ND_all<ND_lim[0]] = np.nan
ND_all[ND_all>ND_lim[1]] = np.nan
F_lim = [-2,1.1]
F_all[F_all<F_lim[0]] = np.nan
F_all[F_all>F_lim[1]] = np.nan

for i, sp in enumerate(sp_names):
    for j in range(1,len(ND_all)):
        if np.isfinite(ND_all[j,i]*F_all[j,i]):
            ax[richness[j]-1].text(ND_all[j,i], F_all[j,i], sp,
                                   color = colors[i],
                               **text_args)
            
for i in range(1,1+max(richness)):
    ax[-1].scatter(ND_all[richness == i], F_all[richness == i], label = "Richness {}".format(i),
                   s = 5*i)

ax[-1].legend()

ax[0].set_xlim(ND_lim)
ax[0].set_ylim(F_lim)
fig.savefig("Figure_ap_assembly_on_NFD.pdf")
