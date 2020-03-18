"""
@author: J.W.Spaak
Create the extended coexstence plot together with bar plots
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches

import plot_annual_plant_definitions as ap



fig, axa = plt.subplots(3,3, figsize = (9,9), sharex = True, sharey = True)
ax = axa.flatten()
plt.xlim(min(ap.interspec), max(ap.interspec))

for i,key in enumerate(ap.keys):
    print(key)
    ND_range = [-0.5,1.5]
    ax[i].plot(ap.interspec, ap.ND[key], label = key, linewidth = 2, alpha = 1, 
                 color = ap.colors[key])
    if i == 0:
        axorg = ax[i].twinx()
        ax2 = axorg
    else:
        ax2 = ax[i].twinx()
        ax2.get_shared_y_axes().join(ax2, axorg)
    ax2.plot(ap.interspec, ap.FD[key], label = key, linewidth = 2, alpha = 1, 
                 color = ap.colors[key], linestyle = "--")
    ax[i].set_title(key)
    rect_facilitation = patches.Rectangle([ap.interspec[0],1],
                -ap.interspec[0], ND_range[1]-1, fill = False, linestyle = ":")
    rect_norm = patches.Rectangle([0,0],ap.sign,1, fill = False)
    rect_comp = patches.Rectangle([ap.sign*1,0],ap.interspec[-1], ND_range[0]
                                  , fill = False, linestyle = "--")
    ax2.set_ylim([-3,3])
    if i%3==2:
        ax2.set_ylabel(r"Fitness differences ($\mathcal{F}$)")
    
    ax[i].add_patch(rect_norm)
    ax[i].add_patch(rect_facilitation)
    ax[i].add_patch(rect_comp)

for i in range(3):
    axa[i,0].set_ylabel(r"Niche differences ($\mathcal{N}$)")
    axa[-1,i].set_xlabel(r"Interspecific interaction ($\alpha$)")
plt.xticks([0,ap.sign*1])
ax[0].set_ylim(*ND_range)

ax[0].set_yticks([0,1])
fig.savefig("ap_annual_plant.pdf")
