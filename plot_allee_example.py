import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

import networkx as nx

from scipy.optimize import minimize
from scipy.integrate import solve_ivp

from NFD_definitions.numerical_NFD import NFD_model, InputError

# number of species
n_spec = 4
species = ["GORO", "LARO", "POAR", "POGN"]

coms = {}

# define species parameters
# all species parameters were chery picked from the bayesian distribution of parameter fits
# they were chosen to show coexistence and intra-specific facilitation
coms["c_matrix"] = np.array([[-0.09900384, -0.00238698, -0.00618794, -0.07544004],
       [-0.05214694, -0.08525456, -0.06602491, -0.06713724],
       [-0.01767611, -0.13315057, -0.25160631, -0.18881568],
       [-0.00302025, -0.14332788, -0.10785323, -0.05064162]])
coms["A_0"] = np.array([[-0.00318617,  0.01072818, -0.00497207,  0.01596994],
        [ 0.09677573, -0.07883874, -0.04840813, -0.00690859],
        [ 0.15020723, -0.07356559, -0.11115578,  0.12129592],
        [ 0.07534394,  0.08312718, -0.14698912, -0.27744169]])
coms["A_slope"] = np.array([[-0.52574907, -0.44019756, -0.12685917, -0.01618502],
       [-0.15818559, -0.10697703, -0.31448185, -0.29972812],
       [-0.44651002, -0.12246011, -0.29592039, -0.19111641],
       [-0.33352455, -0.19490189, -0.71548282, -0.56230952]])
coms["N_opt"] = np.array([[ 3.67685079, 16.93247557,  1.13624294, 17.30170116],
       [ 0.36008801,  0.47877344,  0.80422172,  4.33630166],
       [ 0.56014291,  0.2269895 ,  0.38777543,  1.68272275],
       [ 4.20545803,  0.61476899,  0.12664649,  0.7509636 ]])

# reset parameters
coms["g"] = np.array([0.7       , 0.10408643, 0.1       , 0.10202286])
coms["s"] = np.array([0.16641176, 0.6585376 , 0.68028639, 0.59922572])
coms["lambs"] = np.array([ 91.94858363,  26.26786457, 415.9362581 , 365.3773328 ])

# species order
coms["sp_loc"] = np.arange(n_spec)


def model(S, sp_ind = np.arange(4), ret = "growth"):
    np.exp(S)
    A_slope = coms["A_slope"][sp_ind[:, np.newaxis], sp_ind]
    A_0 = coms["A_0"][sp_ind[:, np.newaxis], sp_ind]
    c_matrix = coms["c_matrix"][sp_ind[:, np.newaxis], sp_ind]
    N_opt = coms["N_opt"][sp_ind[:, np.newaxis], sp_ind]
    lambs = coms["lambs"][sp_ind]
    germination = coms["g"][sp_ind]
    survival = coms["s"][sp_ind]
    N = germination*S # focus on germinated species
    # create respective interaction matrix
    
    A = A_0 + c_matrix*((1 - np.exp(A_slope*(N - N_opt)))
                        /(1 + np.exp(A_slope*(N - N_opt))))
    # remove nan values
    #B = A.ravel()
    #B[np.isnan(B)] = np.ravel(A_0 + c_matrix)[np.isnan(B)]
    if ret != "growth":
        return A
    return (1-germination)*survival + germination*lambs*np.exp(A.dot(N))

# compute niche and fitness differences
pars = {} # dictionary to compute niche and fitness differences, starting estimates
pars["c"] = np.array([[1.00000000e+00, 2.30552544e-04, 8.59300126e-04, 2.18598742e-03],
       [4.33740605e+03, 1.00000000e+00, 1.37126255e+00, 1.15003869e+00],
       [1.16373776e+03, 7.29254949e-01, 1.00000000e+00, 2.67181589e-01],
       [4.57459175e+02, 8.69535961e-01, 3.74277286e+00, 1.00000000e+00]])
# to find equilibrium faster
pars["N_star"] = np.array([[  0.        ,  42.30804828, 129.70334436,  38.78275102],
       [ 44.94285436,   0.        , 279.37710591,   6.48626792],
       [ 43.9788569 , 206.16188413,   0.        , 186.30613119],
       [ 46.26857548, 155.59421515, 208.29198941,   0.        ]])

c_start = 60
pars["c"][0, 2] = c_start
pars = NFD_model(lambda N: np.log(model(N, coms["sp_loc"]))
                 ,  n_spec = 4, pars = pars)

from NFD_definitions.visualize_c_computation import plot_c_computations
fig, ax = plot_c_computations(pars, c_range = 501)
ax[0,0].semilogy()
ax[0,0].set_ylim([0.01, None])
ax[2,0].axvline(c_start)
ax[1,0].grid()
###############################################################################
# plot results

fig = plt.figure(figsize=(9,9))
gs = gridspec.GridSpec(2, 2, figure=fig)

# First three plots
ax = np.empty((2,2), dtype = "object")
ax_NFD = fig.add_subplot(gs[1, 1])
ax_aij = fig.add_subplot(gs[0, 1])
ax_eta = fig.add_subplot(gs[1, 0])

# For ax[1,1], divide it into 2x2 smaller subplots
gs_sub = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0,0])

sub_ax = np.empty((2,2), dtype = "object")
sub_ax[0,0] = fig.add_subplot(gs_sub[0, 0])
sub_ax[0,1] = fig.add_subplot(gs_sub[0, 1]
                              , sharex = sub_ax[0,0], sharey = sub_ax[0,0])
sub_ax[1,0] = fig.add_subplot(gs_sub[1, 0]
                              , sharex = sub_ax[0,0], sharey = sub_ax[0,0])
sub_ax[1,1] = fig.add_subplot(gs_sub[1, 1]
                              , sharex = sub_ax[0,0], sharey = sub_ax[0,0])

#ax_NFD.set_xlim(-0.2, 2)
#ax_NFD.set_ylim([-0.2, 2])
ax_NFD.set_xticks([0,1])
ax_NFD.set_yticks([0,1])
colors = ["cyan", "green", "orange", "purple"]
for i, sp in enumerate(species):
    ax_NFD.plot(pars["ND"][i], pars["F"][i], 'o', label = sp,
                 color = colors[i])
ax_NFD.arrow(pars["ND"][0], 1.2, 0, 0.2, color = colors[0],
              width = 0.01)
ax_NFD.plot(pars["ND"][0], 1.2, 'o', color = colors[0])
ax_NFD.set_ylim([-0.2, 1.5])
ax_NFD.legend()
ax_NFD.set_xlabel("Niche differences", fontsize = 14)
ax_NFD.set_ylabel("Fitness differences", fontsize = 14)
ax_NFD.axhline(0, color = "grey")
ax_NFD.axhline(1, color = "grey")
ax_NFD.axvline(0, color = "grey")
ax_NFD.axvline(1, color = "grey")

ax_NFD.plot(ax_NFD.get_xlim(), ax_NFD.get_xlim(), color = "grey")

###############################################################################
# show species interactions
n_surv = 4 # number of species who can coexist
A_eff = np.full((n_surv,n_surv), np.nan)

for i in range(n_surv):
    # interspecific interaction computed based on invasion growth rate
    A_eff[i] = model(pars["N_star"][i], ret = "matrix")[i]
    
    # intraspecific interaction computed based on no-niche growth rate
    
    A_eff[i,i] = model(n_surv*[np.sum(pars["c"][i]*pars["N_star"][i])],
                        ret = "matrix")[i,i]

# Define species positions in a triangle
positions = {0: (0,0), 1: (0,1), 2: (1, 0), 3: (1,1)}

# Create a directed graph
G_inter = nx.DiGraph()

# Add nodes
for i in range(len(A_eff)):
    G_inter.add_node(i)

weights = np.sign(A_eff)*np.sqrt(np.abs(A_eff))
weights = A_eff
color = np.where(np.sign(A_eff)>0, "r", "b")
# Add edges with weights
for i in range(len(A_eff)):
    for j in range(len(A_eff)):
        if i == j: continue
        G_inter.add_edge(i, j, weight = weights[j,i])

# Define edge colors and widths based on interaction strength
edges = G_inter.edges(data=True)
edge_colors = ["blue" if data["weight"] < 0 else "red" for _, _, data in edges]
weight_multiplier = 10
edge_widths = [weight_multiplier*abs(data["weight"])  for _, _, data in edges]  # Scale for visibility

# Draw the graph
nx.draw(
    G_inter, positions, with_labels=True, node_color="None", 
    edge_color=edge_colors, width=edge_widths, arrows=True, connectionstyle="arc3,rad=0.1",
    ax = ax_aij, labels = {i:species[i] for i in range(4)}
)

# add self-loops
def self_loop(x_center, y_center, linewidth, color, ax):
    radius = 0.1
    ax.annotate("", (x_center + radius, y_center), (x_center, y_center + radius),
                arrowprops=dict(arrowstyle="<|-",
                                shrinkA=10,  # creates a gap between the start point and end point of the arrow
                                shrinkB=0,
                                linewidth=linewidth,
                                connectionstyle="angle,angleB=-90,angleA=180,rad=10",
                                color = color))    
    
    ax.annotate("", (x_center, y_center - radius), (x_center + radius, y_center), 
                arrowprops=dict(arrowstyle="-",
                                shrinkA=0, 
                                shrinkB=0,
                                linewidth=linewidth,
                                connectionstyle="angle,angleB=180,angleA=-90,rad=10",
                                color = color))    
    
    ax.annotate("", (x_center - radius, y_center),  (x_center, y_center - radius), 
                arrowprops=dict(arrowstyle="-",
                                shrinkA=0, 
                                shrinkB=0,
                                linewidth=linewidth,
                                connectionstyle="angle,angleB=-90,angleA=180,rad=10",
                                color = color),
                )    
    ax.annotate("", (x_center, y_center + radius), (x_center - radius, y_center), 
                arrowprops=dict(arrowstyle="-",
                                facecolor="k",
                                linewidth=linewidth,
                                shrinkA=0, 
                                shrinkB=0,
                                connectionstyle="angle,angleB=180,angleA=-90,rad=10",
                                color = color))

dy = 0.2
dpos = [[0, -dy], [0, dy], [0,-dy], [0,dy]]
for i in range(n_surv):
    self_loop(*(np.array(positions[i]) + dpos[i]), weight_multiplier*weights[i,i], color[i,i], ax_aij)
    
ax_aij.set_ylim([-0.3,1.3])

# compute niche and fitness differences of excluded species
c = np.ones((n_spec, n_spec))
c[:n_surv, :n_surv] = pars["c"]
N_mono = np.empty(n_spec) # monoculutre equlibirum densities
for i in range(n_surv):
    pars_sub = NFD_model(lambda N: np.log(model(N, sp_ind = np.array([i,3]))))
    N_mono[i] = pars_sub["N_star"][1,0]
    N_mono[-1] = pars_sub["N_star"][0,1]
    c[i,-1] = pars["c"][0,1]
    c[-1,i] = pars["c"][1,0]

###############################################################################
# plot the intraspecific interaction

log = True
if log:
    N_dens = np.geomspace(5e-1, 2e4, 101)
else:
    N_dens = np.linspace(0, 150, 101)
for i in range(n_spec):
    
    ax_eta.plot(N_dens, [np.log(model(np.insert(np.zeros(3), i, N_i), sp_ind = np.arange(4)))[i]
                          for N_i in N_dens], color = colors[i])
    ax_eta.plot(N_mono[i], 0, 'o', color = colors[i])

    ax_eta.plot(N_dens[0], pars["mu"][i], '^', color = colors[i])
    ax_eta.plot(np.sum(pars["c"][i]*pars["N_star"][i]), pars["eta"][i], '*',
                 color = colors[i])
    print(i, np.sum(pars["c"][i]*pars["N_star"][i]))

ax_eta.plot(np.nan, np.nan, 'ko', label = "Equilibrium density")
ax_eta.plot(np.nan, np.nan, 'k*', label = "No-niche growth rate $\eta_i$")
ax_eta.plot(np.nan, np.nan, 'k^', label = "Intrinsic growth rate $\mu_i$")
ax_eta.legend()
    
ax_eta.set_xlim(N_dens[[0,-1]])
ax_eta.axhline(0, linestyle = "--", color = "k")
ax_eta.set_xlabel("Monoculture density $N_i$", fontsize = 14)
ax_eta.set_ylabel("Monoculture\nGrowth rate", fontsize = 14)


if log:
    ax_eta.semilogx()



for i in range(2):
    sub_ax[i,0].set_ylabel("$A_{ij}(N_j)$", fontsize = 14)
    sub_ax[i,1].yaxis.set_visible(False)
    sub_ax[1,i].set_xlabel("$N_j$", fontsize = 14)
    sub_ax[0,i].xaxis.set_visible(False)
sub_ax[0,0].semilogx()
sub_ax[0,0].set_ylim([-0.4, 0.2])

sub_ax = sub_ax.flatten()
for i in range(n_spec):
    N_eta = np.sum(pars["N_star"][i]*pars["c"][i])
    sub_ax[i].plot(np.sum(pars["N_star"][i]*pars["c"][i]),
                   model(np.insert(np.zeros(3), i, N_eta), ret = "matrix")[i,i]
                   , '.', color = colors[i])
N_comp = np.geomspace(1, sub_ax[0].get_xlim()[1]*2, 101)
for i in range(n_spec):
    for j in range(n_spec):
        sub_ax[i].set_title(species[i])
        sub_ax[i].set_title("ABCD"[i], loc = "left")
        sub_ax[i].plot(N_comp, 
                       [model(np.insert(np.zeros(3), j, N), ret = "matrix")[i,j] for N in N_comp],
                       color = colors[j])
        # add the density of the equilibrium
        sub_ax[i].plot(pars["N_star"][i,j], model(pars["N_star"][i], ret = "matrix")[i,j], '.',
                       color = colors[j])
        
ax_NFD.set_title("G", loc = "left")
ax_NFD.set_title("Coexistence Map")
ax_aij.set_title("E", loc = "left")
ax_aij.set_title("Species interactions")
ax_eta.set_title("F", loc = "left")
ax_eta.set_title("Monoculture growth rates")


fig.tight_layout()

fig.savefig("Figure_alle_example.png")


