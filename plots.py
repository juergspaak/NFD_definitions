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
fig, ax = plt.subplots(1,2, figsize = (10,3.5))
bars = np.array([0.1,0.2,0.3,0.4,0.5,0.4,0.3,0.2,0.1])
x = np.arange(len(bars))

# barplots of dependence on limiting factor
ax[0].bar(x,bars, color = "white", alpha = 1, edgecolor = "black")
ax[0].bar(x+4,bars, color = "black", alpha = 1)
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


ax[1].text(0.5,2/3, r"$\frac{1}{N_i}\frac{dN_i}{dt}=$"
                +r"$\sum_{l=1}^mu_{il}R_l-m_i$",
                        ha= "center", fontsize = 24)

ax[1].text(0.5,1/3, r"$\frac{1}{R_l}\frac{dR_l}{dt}=$"+
            r"$r_l\left(1-\frac{R_l}{K_l}\right)-\sum_{i=1}^nu_{il}N_i$",
                        ha= "center", fontsize = 24)

"""t =("Figure: Dependence of two species \n"
    "(white and black) on limiting factors.\n"
    "In the MacArthur resource model the\n"
    "limiting factors are the resources $R_l$\n"
    "and the dependence is the consuption $u_{il}.$")

from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_weight('bold')
ax[1].text(-0.1,0.5,t,va = "center", fontsize = 14, fontproperties = font)"""
ax[1].axis("off")
fig.savefig("Figure, Limiting factors.pdf")

###############################################################################
# Extended regions
fig = plt.figure(figsize = (7,7))
from matplotlib.text import TextPath

x = np.linspace(-1,3,101)
plt.plot(x,x/(1-x), "black")
plt.axis([x[0],2,-1,4])
plt.ylabel(r"Fitness differences$(-\mathcal{F})$", fontsize = fs)
plt.xlabel(r"Niche differnces $(\mathcal{N})$", fontsize = fs)

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

fig.savefig("Extended Coexistence region.pdf")

###############################################################################
# Experimental setup

fig = plt.figure(figsize = (11,7))
x = y = 0.01
plt.arrow(x,y,0.9,0, head_width = 0.01, color = "black")
plt.arrow(x,y,0,0.9, head_width = 0.01, color = "black")

t1 = 0.7
equi = 0.5

plt.hlines(equi,x,0.9+x, color = "black", linestyles = "dotted")
plt.vlines(t1,y,0.9+x, color = "black", linestyles = "dashed")
plt.axis([0,1,0,1])

fs = 16
plt.ylabel("Density", fontsize = fs)
plt.yticks([equi], [r"$N_i^*$"], fontsize = fs)
plt.xticks([t1], [r"$t_1$"], fontsize = fs)
plt.xlabel(r"Time $t$", fontsize = fs)

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# plot experiment 1
speed = 15
time,dt = np.linspace(x,t1,10, retstep = True)
growth = equi/(1+np.exp(-speed*time)*(equi/x-1))
plt.plot(time, growth, 'o', color = "black")
plt.text(0.3,0.33,"1. Monoculture growth\nfrom low abundance",
         va = "top", ha = "center", fontsize = 16)

decline = equi/(1+np.exp(-speed*time)*(equi/0.9-1))+x
plt.scatter(time, decline,facecolor = "none", color = "black")
plt.text(0.3,0.7,"2. Monoculture growth\nfrom high abundance",
         va = "top", ha = "center", fontsize = 16)

time2 = t1+np.arange(3)*dt
plt.plot(time2,equi*np.ones(3),'v', color = "black")
invasion = equi/(1+np.exp(-speed*(time2-t1))*(equi/x-1))+x
plt.plot(time2, invasion, '^', color = "black")

plt.text(0.82,0.7,"3. Invasion\nexperiment;\nresident density",
         va = "top", ha = "center", fontsize = 16)

plt.text(0.82,0.3,"3. Invasion\nexperiment;\ninvader density",
         va = "top", ha = "center", fontsize = 16)

fig.savefig("Experimental_cartoon.pdf")