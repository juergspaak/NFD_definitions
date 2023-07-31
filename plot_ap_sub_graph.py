import numpy as np
import matplotlib.pyplot as plt



from nfd_definitions.numerical_NFD import NFD_model
import pandas as pd

import scheme_plot as schp
    
###############################################################################
# 
A = pd.read_csv("Geijzendorffer 2011.csv", index_col = 0).values
n = len(A)
sp_names = np.arange(len(A)) + 1
# remove species 3 from the community
no_3 = np.where(np.arange(len(A)) != 2)[0]
A = A[no_3[:,np.newaxis], no_3]
sp_names = sp_names[no_3]



fig, ax = plt.subplots(2,2, figsize = (9,9))

ax[0,0].set_title("Focus on species 4")
ax[0,1].set_title("Focus on species 6")

ind = np.where(sp_names != 4)[0]
pars_4 = schp.do_all_LV(A[ind[:,np.newaxis], ind], axc = ax[0,0], plot_non_stable=True,
               sp_names = sp_names[ind])[-1]

ind = np.where(sp_names != 6)[0]
pars_6 = schp.do_all_LV(A[ind[:,np.newaxis], ind], axc = ax[0,1], plot_non_stable=True,
               sp_names = sp_names[ind])[-1]

inv_scheme, equi_all = schp.compute_scheme_for_LV(A)
equi_all = equi_all[np.all(equi_all>-1e-10, axis = 1)]
pars_all = schp.compute_NFD_LV(A)
ax[1,0].set_title("Full community")
a = schp.plot_scheme(inv_scheme, ax[1,0], plot_non_stable = False,
                     sp_names = sp_names)


text_args = dict(ha = "center", va = "center", fontsize = 16)

r_i = inv_scheme
mu_i = np.ones(r_i.shape)
c_ij = pars_all["c"]
eta_i = mu_i - np.diag(A)*np.einsum("ij,nj->ni", c_ij, equi_all)

ND_all = (r_i - eta_i)/(mu_i - eta_i)
F_all = -eta_i/(mu_i - eta_i)

id_1246 = -4
id_126 = -9

id_125 = -10

sp = 5
NFD = np.array([ND_all[[id_126,id_1246],sp_names == sp],
                F_all[[id_126,id_1246],sp_names == sp]]).T
ax[1,1].text(*NFD[0],
             sp, **text_args)
ax[1,1].text(*NFD[1],
             sp, **text_args)
ax[1,1].arrow(*(NFD[0] + 0.05*(NFD[1]-NFD[0])), *(0.9*(NFD[1]-NFD[0])), 
              length_includes_head = True,
              width = 0.01, color = "k")

sp = 4
NFD = np.array([ND_all[[id_125,id_126],sp_names == sp],
                F_all[[id_125,id_126],sp_names == sp]]).T
ax[1,1].text(*NFD[0],
             sp, **text_args)
ax[1,1].text(*NFD[1],
             sp, **text_args)
ax[1,1].arrow(*(NFD[0] + 0.1*(NFD[1]-NFD[0])), *(0.8*(NFD[1]-NFD[0])), 
              length_includes_head = True,
              width = 0.01, color = "k")

for a in ax[[0,0,1],[0,1,0]]:
    a.set_xticks([])
    a.set_ylabel("Species richness")

ax[1,1].set_xlabel("Niche differences")
ax[1,1].set_ylabel("Fitness differences")



ax[1,1].plot([0,1],[0,1], 'r')
ax[1,1].set_xlim([0,2.5])
ax[1,1].axvline(0, color = "grey")
ax[1,1].axvline(1, color = "grey")
ax[1,1].axhline(0, color = "grey")
ax[1,1].axhline(1, color = "grey")

ax[1,1].set_title("Effect of focal species")

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCD"[i], loc = "left")

fig.savefig("Figure_ap_subgraph.pdf")
