"""
@author: J.W.Spaak
Create the plots used in the manuscript
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'Times New Roman'
###############################################################################
# Example of Mac Arthur resource model with NO
fig, ax = plt.subplots(1,2, figsize = (10,3.5), sharey = True)
bars = np.array([0.1,0.2,0.3,0.4,0.5,0.4,0.3,0.2,0.1])
x = np.arange(len(bars))

# barplots of dependence on limiting factor
ax[0].set_title("A")
ax[0].bar(x,bars, color = "white", alpha = 1, edgecolor = "black")
ax[0].bar(x+4,bars, color = "grey", alpha = 1, edgecolor = "black")
ax[0].bar(x,bars, color = "white", alpha = 0.5, edgecolor = "black")

fs = 14
#ax[0].text(6.4,0.04,"Niche\nOverlap", ha = "center",fontsize = fs-2,
#        bbox = {"facecolor": "blueviolet","boxstyle": "round"})

ax[0].set_xticks([])
ax[0].set_yticks([])

ax[0].set_xlabel("Resource $R_l$\n(Limiting factor)", fontsize = fs)
ax[0].set_ylabel("Consumption $u_{il}$\n(Dependence on limiting factor)"
           , fontsize = fs)

ax[0].axis([-1,13,None,0.6])

# multispecies community
# barplots of dependence on limiting factor
ax[1].set_title("B")
ax[1].bar(x,bars, color = "white", alpha = 1, edgecolor = "black")
ax[1].bar(x+12,bars, color = "black", alpha = 1, edgecolor = "black")
ax[1].bar(x+4,bars, color = "lightgrey", alpha = 1, edgecolor = "black")
ax[1].bar(x,bars, color = "white", alpha = 0.5, edgecolor = "black")
ax[1].bar(x+12,bars, color = "black", alpha = 0.5, edgecolor = "black")

ax[1].set_xlabel("Resource $R_l$\n(Limiting factor)", fontsize = fs)
ax[1].set_xticks([])
fig.savefig("Figure_Limiting_factors.eps")

###############################################################################
# Extended regions
fig = plt.figure(figsize = (7,7))
from matplotlib.text import TextPath

x = np.linspace(-1,3,101)
plt.plot(x,x/(1-x), "black")
plt.axis([x[0],2,-1,4])
plt.ylabel(r"Fitness differences$(-\mathcal{F})$", fontsize = fs)
plt.xlabel(r"Niche differences $(\mathcal{N})$", fontsize = fs)

plt.axhline(y=0, color = "black", linestyle = ":")
plt.axvline(x=0, color = "black", linestyle = ":")
plt.axvline(x=1, color = "black", linestyle = ":")

plt.xticks([0,1])
plt.yticks([0,-1])

ms = 12
offset = (-4,-4)
# plot the varius definitions
plt.plot([-0.5,-0.5], [-0.1,0.1], linestyle = "",  markersize = ms,
         marker = TextPath(offset, "1"), color = "black",
         label = "priority effects")

plt.plot([0,0], [0,0], '>',  linestyle = "",  markersize = ms,
         marker = TextPath(offset, "2"), color = "black", label = "neutrality")

plt.plot([0.232,0.232], [-0.12968,0.06468], linestyle = "",  markersize = ms,
         marker = TextPath(offset, "3"),
         color = "black", label = "coexistence\n(see exp. Fig. 4)") 

plt.plot([0.5,0.5], [-0.3,3], linestyle = "",  markersize = ms,
         marker = TextPath(offset, "4"), color = "black",
         label = "competitive\nexclusion")

plt.plot([1.2,0.8], [2,-0.8], 's', linestyle = "",  markersize = ms,
         marker = TextPath(offset, "5"), color = "black", label = "parasitism")

plt.plot([1.4,1.4], [-0.9,1.7], linestyle = "",  markersize = ms,
         marker = TextPath(offset, "6"), color = "black", label = "mutualism")

coex_x = np.append(np.linspace(x[0],1-1e-3),x[[-1,-1,0]])
coex_y = coex_x/(1-coex_x)
coex_y[-1] = -1
coex_y[np.isinf(coex_y)] = 10**10
plt.fill(coex_x,coex_y, color = "grey",alpha = 0.5)

plt.legend(numpoints = 1)

fig.savefig("Extended_Coexistence_region.eps")