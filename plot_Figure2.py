"""
@author: J.W.Spaak
Create Fiugre 2
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'Times New Roman'
###############################################################################
# Example of Mac Arthur resource model with NO
fig, ax = plt.subplots(1,3, figsize = (10,3.5), sharey = True)
bars = np.array([0.1,0.2,0.3,0.4,0.5,0.4,0.3,0.2,0.1])
x = np.arange(len(bars)+12)

# barplots of dependence on limiting factor, unscaled
ax[0].set_title("A")
sp = np.zeros((2,len(x)))
sp[0,:len(bars)] = bars
sp[1,4:(4+len(bars))] = 0.6*bars

ax[0].bar(x,sp[0], color = "white", alpha = 1, edgecolor = "black")
ax[0].bar(x,sp[1], color = "grey", alpha = 1, edgecolor = "black")
ax[0].bar(x, np.amin(sp, axis = 0), color ="lightgrey", edgecolor = "black")

ax[0].set_xticks([])
ax[0].set_yticks([])


sp = np.zeros((2,len(x)))
sp[0,:len(bars)] = bars
sp[1,4:(4+len(bars))] = bars
# barplots of dependence on limiting factor
ax[1].set_title("B")
ax[1].bar(x,sp[0], color = "white", alpha = 1, edgecolor = "black")
ax[1].bar(x,sp[1], color = "grey", alpha = 1, edgecolor = "black")
ax[1].bar(x, np.amin(sp, axis = 0), color ="lightgrey", edgecolor = "black")

fs = 14

ax[1].set_xticks([])
ax[1].set_yticks([])

ax[0].set_xlabel("Resource $R_l$\n(Limiting factor)", fontsize = fs)
ax[1].set_xlabel("Resource $R_l$\n(Limiting factor)", fontsize = fs)
ax[0].set_ylabel("Consumption $u_{il}$\n(Dependence on limiting factor)"
           , fontsize = fs)


ax[0].axis([-1,13,None,0.6])
ax[1].axis([-1,13,None,0.6])

# multispecies community
# barplots of dependence on limiting factor

sp = np.zeros((3,len(x)))
sp[0,:len(bars)] = bars
sp[1,4:(4+len(bars))] = bars
sp[2,12:(12+len(bars))] = bars

ax[2].set_title("C")
ax[2].bar(x,sp[0], color = "white", edgecolor = "black",
  linewidth = 0.5)
ax[2].bar(x,sp[2], color = "black", edgecolor = "black",
  linewidth = 0.5)
ax[2].bar(x,sp[1], color = "grey", edgecolor = "black",
  linewidth = 0.5)
ax[2].bar(x,np.amin(sp[:2], axis = 0), color = "lightgrey", edgecolor = "black",
  linewidth = 0.5)
ax[2].bar(x,np.amin(sp[1:], axis = 0), color = (0.25, 0.25, 0.25)
  , edgecolor = "black",
  linewidth = 0.5)

ax[2].set_xlabel("Resource $R_l$\n(Limiting factor)", fontsize = fs)
ax[2].set_xticks([])
fig.savefig("Figure2.pdf")