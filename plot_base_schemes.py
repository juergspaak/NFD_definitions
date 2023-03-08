import numpy as np
import matplotlib.pyplot as plt

import scheme_plot as schp

fig, ax = plt.subplots(2,2, sharex = True, sharey = True,
                       figsize = (9,9))
ax = ax.flatten()

n = 3
sp_names = np.arange(1, n+1)
# full community
A = np.eye(3)
r_i, equi_all, inv_graph, acyclic, permanent = schp.do_all(A, axc = ax[0], sp_names = sp_names)

# remove entry with species 1 and 1 coexisting
r_i = r_i[list(range(5)) + [-2,-1]]
schp.plot_scheme(r_i, axc = ax[1], sp_names = sp_names)

A = np.array([[1,0,0],
              [0,1,0],
              [2,2,1]])
r_i, equi_all, inv_graph, acyclic, permanent = schp.do_all(A, axc = ax[2], sp_names = sp_names)
A = np.array([[1,0,2],
              [2,1,0],
              [0,2,1]])
r_i, equi_all = schp.compute_scheme_for_LV(A)
r_i = r_i[:-1]
schp.plot_scheme(r_i, axc = ax[3], sp_names = sp_names)


###############################################################################
# add layout
ax[0].set_ylabel("Species richness")
ax[2].set_ylabel("Species richness")

ax[0].set_title("Traditional community")
ax[1].set_title("Not all N-1 form\nbut all coexist")
ax[2].set_title("Not all coexist\nwe analyse a subcommunity")
ax[3].set_title("Cyclic subgraph")

ax[0].set_xticks([])

for i, a in enumerate(ax):
    a.set_title("ABCD"[i], loc = "left")
    
fig.savefig("Figure_base_schemes.pdf")