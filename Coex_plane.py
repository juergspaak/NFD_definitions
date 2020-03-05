"""
@author: J.W.Spaak
Create the extended coexstence plot together with bar plots
"""

import matplotlib.pyplot as plt
import numpy as np
from string import ascii_lowercase as letters

import matplotlib as mpl
from matplotlib.textpath import TextPath
from numerical_NFD import NFD_model

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
plt.rcParams["font.family"] = 'Times New Roman'

def LV_model(N,A):
    # lotka volterra model
    return 1 - A.dot(N)

# interaction matrices
interaction_matrices = np.array([[[ 1.0, 1.2], [ 1.4, 1.0]],
                                 [[ 1.0, 1.0], [ 1.0, 1.0]],
                                 [[ 1.0, 0.3], [ 1.2, 1.0]],
                                 [[ 1.0, 0.4], [-0.2, 1.0]],
                                 [[ 1.0,-0.6], [-0.4, 1.0]]])
    
c_help = 0.5
A = np.array([[0.25*c_help, 0.3], [0.3125*c_help, 1.5]])/2
interaction_matrices[2] = A
    
ND, FD, f0, fc, ri, c = np.empty((6, len(interaction_matrices), 2))

for i, A in enumerate(interaction_matrices):
    pars = NFD_model(LV_model, args = (A,))
    ND[i] = pars["ND"]
    FD[i] = pars["FD"]
    f0[i] = pars["f0"]
    fc[i] = pars["fc"]
    ri[i] = pars["r_i"]
    c[i] = pars["c"][[1,0],[0,1]]

###############################################################################
# Extended regions
    
fig, ax = plt.subplots(2,2, sharey = "row", sharex = "row",
                       figsize = (8,8))

viridis = mpl.cm.get_cmap('viridis')
colors = viridis(np.linspace(0,1,len(interaction_matrices)))
# use cathegorical colors from https://jfly.uni-koeln.de/html/color_blind/
colors = np.array([[0,0,0], #black
          [230,159,0], # orange
          [86,180,233], # sky blue
          [0,158,115], # bluish green
          [240,228,66], # yellow
          [0,114,178],# blue
          [213,94,0], # vermillion
          [204,121,167]])/255 # reddish purple
colors = colors[1:]
labels = ["priority effects", "neutrality", "competitive\nexclusion",
          "parasitism", "mutualism"]
offset = (-4,-4)
dist = np.arange(2)*(len(interaction_matrices)+1)
lw = 0.5
xticks = []
for i in range(2):
    # plot the extended coexistence range
    
    # background
    x = np.linspace(-1,2,100)
    ax[0,i].plot(x,x/(1-x), "black")
    ax[0,i].set_ylim([-1,2])
    ax[0,i].set_xlim([-1,2])
    ax[0,i].axhline(y=0, color = "black", linestyle = ":")
    ax[0,i].axvline(x=0, color = "black", linestyle = ":")
    ax[0,i].axvline(x=1, color = "black", linestyle = ":")
    
    coex_x = np.append(np.linspace(x[0],1-1e-3),x[[-1,-1,0]])
    coex_y = coex_x/(1-coex_x)
    coex_y[-1] = -1
    coex_y[np.isinf(coex_y)] = 10**10
    ax[0,i].fill(coex_x,coex_y, color = "grey",alpha = 0.5)
    
    # ND and FD values
    for j in range(len(interaction_matrices)):
        ax[0,i].plot(ND[j,i], -FD[j,i], marker = TextPath(offset, letters[j]),
          color = colors[j], markersize = 12, label = labels[j],
          linestyle = "")
        
    ax[1,i].axhline(y = 1, color = "black", linestyle = "-", linewidth = lw)
    ax[1,i].axhline(y = 0, color = "black", linestyle = "-", linewidth = lw)
    
    for j in range(len(interaction_matrices)):
        ax[1,i].bar(dist + j, [ri[j,i], fc[j,i]],
          color = colors[j])
        xticks.extend(dist+j)
    ax[1,i].set_xlim([-1,-1 + 2*(len(interaction_matrices) +1)])
    ax[1,i].axvline(x = len(interaction_matrices),
      color = "black", linestyle = "-", linewidth = 1.5)
    ax[1,i].axvline(x = 2*len(interaction_matrices) +1,
      color = "black", linestyle = "-", linewidth = 1.5)
    ax[1,i].text(ax[1,i].get_xlim()[0]+0.1, 1, "intrinsic\ngrowth rate",
      va = "bottom", ha = "left", fontsize = 12)

# add bars for neutrality for visibility
j = 1
val = 0.01
ax[1,0].bar(dist + j, 2*[val],
          color = colors[j])
ax[1,0].bar(dist + j, 2*[-val],
          color = colors[j])
ax[1,1].bar(dist + j, 2*[val],
          color = colors[j])
ax[1,1].bar(dist + j, 2*[-val],
          color = colors[j])
       
# add data from experiment
try:
    pars_exp
except NameError:
    from Exp_plots import pars as pars_exp
ax[0,0].plot(pars_exp["ND"][0], -pars_exp["FD"][0], 'k^',
  label = "Cyanobacteria\nexperiment")
ax[0,1].plot(pars_exp["ND"][1], -pars_exp["FD"][1], 'k^')
ax[0,0].legend(fontsize = 10, loc = "upper left")

# add axis layouts
fs = 16
ax[0,0].set_ylabel(r"Fitness differences$(-\mathcal{F})$", fontsize = fs)
ax[0,0].set_xlabel(r"Niche differences$(\mathcal{N})$", fontsize = fs)
ax[0,1].set_xlabel(r"Niche differences$(\mathcal{N})$", fontsize = fs)

# add ticks
ax[0,0].set_xticks([0,1])
ax[0,0].set_yticks([-1,0])

# growth rate panels
xlim = ax[1,0].get_xlim()
for i in range(2):
    ax_copy = ax[1,i].twiny()
    ax_copy.set_xticks(sorted(set(xticks)))
    ax_copy.set_xlim(xlim)
    ax_copy.set_xticklabels(2*letters[:len(interaction_matrices)])


ax[1,0].set_xticks([xlim[0]*3/4 + xlim[1]*1/4,
                   xlim[0]*1/4 + xlim[1]*3/4])
ax[1,0].set_xticklabels(["invasion\ngrowth rate", "no-niche\ngrowth rate"]
    , fontsize = 12)
ax[1,1].set_xticklabels(["invasion\ngrowth rate", "no-niche\ngrowth rate"]
    , fontsize = 12)
    
ax[1,0].set_yticks([0,1])
ax[1,0].set_ylabel("per capita\ngrowth rate", fontsize = fs)

# panel titles
ax[0,0].set_title("A: Species i", loc = "left")
ax[0,1].set_title("B: Species j", loc = "left")
ax[1,0].set_title("C: Species i", loc = "left")
ax[1,1].set_title("D: Species j", loc = "left")
fig.tight_layout()

fig.savefig("Extended_Coexistence_region.pdf")
fig.savefig("Extended_Coexistence_region.eps")