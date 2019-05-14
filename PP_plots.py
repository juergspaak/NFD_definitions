"""
@author: J.W.Spaak
Create plots of manusrcipt for PP
"""

import matplotlib.pyplot as plt
import numpy as np

import matplotlib as mpl
label_size = 18
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size 

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
plt.style.use('dark_background')

###############################################################################
# Example of Mac Arthur resource model with NO
def save_res(counter = [1]):
    print(counter[0])
    plt.savefig("PP_slides/Resources_example_{}.png".format(counter[0]),
                transparent = "True")
    counter[0] += 1

fig, ax = plt.subplots(2,2, figsize = (9,7), sharex = True, sharey = True)
bars = np.array([0.1,0.2,0.3,0.4,0.5,0.4,0.3,0.2,0.1])
x = np.arange(len(bars))

ax[0,0].set_xticks([])
ax[0,0].set_yticks([])

ax_all = fig.add_subplot(1,1,1, frameon = False)
ax_all.set_xticks([])
ax_all.set_yticks([])

fs = 22
ax_all.set_xlabel("Resource $R_l$\n(Limiting factor)", fontsize = fs)
ax_all.set_ylabel("Consumption $u_{il}$\n(Dependence on limiting factor)"
           , fontsize = fs)

axs = ax.flatten()
dist = [4,4,10,0]
size = [1,0.5,1,1]
NO = [r"$\mathcal{N} =\frac{2}{3}$", r"$\mathcal{N}_r = \frac{3}{4}$" + "\n"
   r"$\mathcal{N}_b = \frac{1}{2}$", r"$\mathcal{N} = 0$", r"$\mathcal{N} = 1$"]
ax[0,0].set_xlim(-7,15)
ax[0,0].set_ylim(0,0.6)
for i in range(4):
    # barplots of dependence on limiting factor
    axs[i].bar(x - dist[i]/2,bars, color = "red", alpha = 1)
    axs[i].bar(x + dist[i]/2,size[i]*bars, color = "blue", alpha = 1)
    axs[i].bar(x - dist[i]/2,bars, color = "red", alpha = 0.5)
    axs[i].text(-7, 0.55, NO[i], fontsize = fs, verticalalignment = "top")
    save_res()


###############################################################################
def save_coex(counter = [1]):
    print(counter[0])
    plt.savefig("PP_slides/Coexistence region_{}.png".format(counter[0]),
                transparent = "True")
    counter[0] += 1    
    
# Extended regions
fig = plt.figure(figsize = (9,9))

x = np.linspace(-1,3,101)
plt.plot(x,x/(1-x), "white", linewidth = 4)
plt.axis([x[0],2,-1,4])
plt.ylabel(r"Fitness differences$(-\mathcal{F})$", fontsize = fs)
plt.xlabel(r"Niche differnces $(\mathcal{N})$", fontsize = fs)

plt.axhline(y=0, color = "white", linestyle = ":", linewidth = 4)
plt.axvline(x=0, color = "white", linestyle = ":", linewidth = 4)
plt.axvline(x=1, color = "white", linestyle = ":", linewidth = 4)

plt.xticks([0,1])
plt.yticks([0,-1])


coex_x = np.append(np.linspace(x[0],1-1e-3),x[[-1,-1,0]])
coex_y = coex_x/(1-coex_x)
coex_y[-1] = -1
coex_y[np.isinf(coex_y)] = 10**10

save_coex()
plt.fill(coex_x,coex_y, color = "white",alpha = 0.5)
save_coex()
ms = 10
# plot the varius definitions
plt.plot([-0.5,-0.5], [-0.1,0.1], 'p',  markersize = ms,
         color = "white", label = "priority effects")

plt.plot([0,0], [0,0], '>',  markersize = ms,
         color = "white", label = "neutrality")

plt.plot([0.2,0.2], [-0.3,3], 'D', markersize = ms,
         color = "white", label = "competitive\nexclusion")

plt.plot([0.6,0.6], [-0.5,1.2], '*',  markersize = ms,
         color = "white", label = "stable\ncoexistence")
plt.legend(numpoints = 1, fancybox = True, framealpha = 0.5,
           fontsize = 16, loc = "upper left")
save_coex()

plt.plot([1.2,0.8], [2,-0.8], 's', markersize = ms,
         color = "white", label = "parasitism")

plt.plot([1.4,1.4], [-0.9,1.7], 'P',  markersize = ms,
         color = "white", label = "mutualism")
plt.legend(numpoints = 1, fancybox = True, framealpha = 0.5, 
           loc = "upper left", fontsize = 16)
save_coex()


