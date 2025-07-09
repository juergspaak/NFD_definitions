"""
@author: J.W. Spaak, j.w.spaak@gmail.com

"""
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


matplotlib.rcParams.update({'font.size': 22})

fig = plt.figure(figsize = (8,8))
ax = plt.gca()
plt.axis([-1,2,-1,2])
ax.set_xticks([0,1])
ax.set_yticks([0,1])

ax.axhline(0, color = "k")
ax.axhline(1, color = "k")
ax.axvline(0, color = "k")
ax.axvline(1, color = "k")

ax.plot([-1,2], [-1,2], "k")

fs_small = 13
col_line = "black"
eq = "grey"
col_area = "darkblue"
###############################################################################
# labels for FD:
ax.text(-1,1.5, "Negative\nequilibrium\nabundance\n", 
        horizontalalignment = "center", fontsize = fs_small, color = col_area,
        rotation = 90, verticalalignment = "center",
        multialignment = "center")
ax.text(-1,0.5, "Lower equilibrium\nabundance\nthan competitor\n", 
        horizontalalignment = "center", fontsize = fs_small, color = col_area,
        rotation = 90, verticalalignment = "center",
        multialignment = "center")

ax.text(-1,-0.5, "Higher equilibrium\nabundance\nthan competitor\n", 
        horizontalalignment = "center", fontsize = fs_small, color = col_area,
        rotation = 90, verticalalignment = "center",
        multialignment = "center")


###############################################################################
# labels for division lines

x_axis_label = 1.5
ax.text(x_axis_label,0,"$\eta_i=0$",
        verticalalignment = "bottom", horizontalalignment = "center",
        fontsize = fs_small, color = eq)
ax.text(x_axis_label,0,"$N_i^*=c_{ij}N_j^*$",
        verticalalignment = "top", horizontalalignment = "center",
        fontsize = fs_small, color = col_line)

ax.text(x_axis_label,1,"$\mu_i=0$",
        verticalalignment = "bottom", horizontalalignment = "center",
        fontsize = fs_small, color = eq)
ax.text(x_axis_label,1,"$N_i^*=0$",
        verticalalignment = "top", horizontalalignment = "center",
        fontsize = fs_small, color = col_line)

y_axis_labels = 1.5
ax.text(1,y_axis_labels,"No species\ninteractions",
        verticalalignment = "center"
        ,horizontalalignment = "right", fontsize = fs_small,
        color = col_line, rotation = 90)
ax.text(1,y_axis_labels,"$r_i=\mu_i$",
        verticalalignment = "center"
        ,horizontalalignment = "left", fontsize = fs_small,
        color = eq, rotation = 90)

ax.text(0,y_axis_labels,"No frequency\ndependence",
        verticalalignment = "center"
        ,horizontalalignment = "right", fontsize = fs_small,
        color = col_line, rotation = 90)
ax.text(0,y_axis_labels,"$r_i=\eta_i$",
        verticalalignment = "center"
        ,horizontalalignment = "left", fontsize = fs_small,
        color = eq, rotation = 90)


ax.text(1.6,1.6, "$r_i=0$\n", rotation = 45,
        horizontalalignment = "center",
        verticalalignment = "center", fontsize = fs_small+2, color = eq)

###############################################################################
# labels for ND:

ax.text(-0.5,-1,"\nStronger inter-\nthan intraspecific\ninteractions", verticalalignment = "center"
        ,horizontalalignment = "center", fontsize = fs_small,
        color = col_area, rotation = 0)
ax.text(0.5,-1,"\nWeaker inter-\nthan interspecific\ninteractions", verticalalignment = "center"
        ,horizontalalignment = "center", fontsize = fs_small,
        color = col_area, rotation = 0)
ax.text(1.5,-1,"\nInterspecific\nand intrapspecific\ndiffer in sign", verticalalignment = "center"
        ,horizontalalignment = "center", fontsize = fs_small,
        color = col_area, rotation = 0)



ax.text(0.5,0.5+0.02, "Persistence\nline", rotation = 45,
        horizontalalignment = "center",
        verticalalignment = "center", fontsize = fs_small, color = col_line)

ax.set_xlabel(r"$\mathcal{N}_i=\frac{r_i-\eta_i}{\mu_i-\eta_i}$")
ax.set_ylabel(r"$\mathcal{F}_i=\frac{0-\eta_i}{\mu_i-\eta_i}$")

# label each panel with the corresponding equation

def labels(loc, order):
    ax.text(*loc, "${} \gtreqless {}$".format(*growth_rates[order[:2]]),
            ha = "center", va = "bottom", fontsize = fs_small+2,
            color = "r")
    ax.text(*loc, "$\gtreqless {}\gtreqless {}$".format(*growth_rates[order[2:]]),
            ha = "center", va = "top", fontsize = fs_small+2, color = "r")

growth_rates = np.array(["\mu_i", "r_i", "\eta_i", "0"])
labels([0.5, -0.5], [0,1,2,3])
labels([1.5, -0.5], [1,0,2,3])
labels([-0.25, -0.75], [0,2,1,3])
labels([-0.75, -0.25], [0,2,3,1])

labels([0.75, 0.25], [0,1,3,2])
labels([0.25, 0.75], [0,3,1,2])
labels([-0.5, 0.5], [0,3,2,1])
labels([1.5, 0.5], [1,0,3,2])

labels([0.5, 1.5], [3,0,1,2])
labels([-0.5, 1.5], [3,0,2,1])
labels([1.25, 1.75], [3,1,0,2])
labels([1.75, 1.25], [1,3,0,2])


fig.savefig("Figure_allee_NFD_map_explained.pdf")

matplotlib.rcParams.update({'font.size': 10.0})