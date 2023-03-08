import numpy as np
import matplotlib.pyplot as plt

from numerical_NFD import NFD_model
import pandas as pd

import scheme_plot as schp

fig, ax = plt.subplots(2,3, figsize = (9,9))

A = pd.read_csv("Vandermeer_1969.csv", index_col = 0).values
n = len(A)

# coexisting subcommunity
coex = np.arange(2)

# plot invasion scheme
schp.do_all_LV(A, axc = ax[0,0], sp_names = "index_at_1",
            plot_min_i_com = False)
schp.do_all_LV(A[coex, coex[:,np.newaxis]], axc = ax[0,1], sp_names = coex+1,
            plot_min_i_com = False)

equi = np.zeros(n)
equi[coex] = np.linalg.solve(A[coex[:,np.newaxis], coex], np.ones(len(coex)))

# compute niche and fitness differences
pars = NFD_model(lambda N: 1 - A[coex[:,np.newaxis], coex].dot(N), len(coex))


ND_old = 1 - np.sign(A)*np.sqrt(np.abs(A*A.T))
F_old = 1 - np.sqrt(np.abs(A.T/A))
# compute c for all species
c = np.sqrt(np.abs(A/A.T))
# compute niche and fitness differences for competitively exclused species
mu_i = 1
r_i = 1 - A.dot(equi)
eta_i = 1 - c.dot(equi)
N_new = (r_i - eta_i)/(mu_i-eta_i)
F_new = -eta_i/(mu_i-eta_i)

other = np.full(n, True, dtype = bool)
other[coex] = False
ax[0,-1].plot(pars["ND"], pars["F"], 'r^')
ax[0,-1].plot(N_new[other], F_new[other], 'bo')

com = [0,3]
ax[0,-1].plot(ND_old[com, com[::-1]], F_old[com, com[::-1]], '^', color = "purple")
com = [2,3]
ax[0,-1].plot(ND_old[com, com[::-1]], F_old[com, com[::-1]], '^', color = "orange")
    
###############################################################################
# 
A = pd.read_csv("Geijzendorffer 2011.csv", index_col = 0).values
A[A==0] = 0.01 #♣ to avoid division by zero
n = len(A)

pars = schp.compute_NFD_LV(A)
case = pars["ND"]<pars["F"]

ax[1,-1].plot(pars["ND"][case], pars["F"][case], 'bo')
ax[1,-1].plot(pars["ND"][~case], pars["F"][~case], 'ro')

scheme, equis = schp.compute_scheme_for_LV(A)
schp.plot_scheme(scheme, ax[1,0], plot_non_stable=False, fs = 7,
                 sp_names = "index_at_1", plot_min_i_com=False)

pres = np.where(~case)[0]
A_foc = A[pres[:,np.newaxis], pres]
scheme, equis = schp.compute_scheme_for_LV(A_foc)
schp.plot_scheme(scheme, ax[1,1], plot_non_stable=False, fs = 7,
                 sp_names = np.array([1,2,4,5,6]), plot_min_i_com=False)

# store matrix for Sebastian schreiber
import pandas as pd
df = pd.DataFrame(A, index = ["sp {}".format(i) for i in range(len(A))],
                  columns =  ["sp {}".format(i) for i in range(len(A))])
df.to_csv("Geijzendorffer 2011.csv")


###############################################################################
# add layout

for i in range(2):
    ax[i,-1].axvline(0, color = 'k')
    ax[i,-1].axvline(1, color = 'k')
    ax[i,-1].axhline(0, color = 'k')
    ax[i,-1].axhline(1, color = "k")
    ax[i,-1].plot([0,1], [0,1], 'r--')
    
for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEF"[i], loc = "left")
    if i%3 != 2:
        a.set_ylabel("Species richness")
        a.set_xticks([])

ax[0,0].set_title("Full community\ninvasion graph")
ax[0,1].set_title("Permanent subcommunity\ninvasion graph")
ax[0,2].set_title("Explaining\ncoexistence")


ax[0,0].set_xticks([])
ax[0,1].set_xticks([])

"""ax[0,0].set_yticks(np.arange(5))
ax[0,1].set_yticks(np.arange(3))
ax[1,0].set_yticks(np.arange(7))
ax[1,1].set_yticks(np.arange(6))"""

ax[0,2].set_ylabel("Fitness differences")
ax[1,2].set_ylabel("Fitness differences")
ax[1,2].set_xlabel("Niche differences")


fig.tight_layout()
    
fig.savefig("Figure_NFD_example.pdf")