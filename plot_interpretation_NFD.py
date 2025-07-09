import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2,2, figsize = (9,9), sharex = False, sharey = False)

mu = 1
eta = -0.4
r = 0.6

N = (r-eta)/(mu - eta)

ylim = np.array([-0.1, 1.1])
xlim = np.array([-0.5, 1.1])

# interpolations
ax[0,0].plot([eta, eta, xlim[0]], [ylim[0], 0, 0], "k")
ax[0,0].plot([mu, mu, xlim[0]], [ylim[0], 1, 1], 'k')
ax[0,0].plot([r, r, xlim[0]], [ylim[0], N, N], 'k--')

ax[0,0].plot([eta, mu], [0, 1], 'r')

ax[0,0].set_ylim(ylim)
ax[0,0].set_xlim(xlim)

ax[0,0].set_xticks([eta, r, mu])
ax[0,0].set_xticklabels(["$\eta_i$", "$r_i$", "$\mu_i$"], fontsize = 16)

ax[0,0].set_yticks([0, N, 1])
ax[0,0].set_yticklabels([0, r"$\frac{r_i-\eta_i}{\mu_i-\eta_i}$", 1], fontsize = 16)

ax[0,0].set_xlabel("(hypothetical)\ninvasion growth rates", fontsize = 14)
ax[0,0].set_ylabel("Niche differences\n$\mathcal{N}_i$", fontsize = 14)
ax[0,0].set_title("Definition of\nNiche differences")

###############################################################################
# interpretation of niche differences

colors = ["darkblue", "blue", "green"]

xlim = np.array([-1.8,2.4])
ylim = np.array([-1,2])

# interpolations
ax[0,1].plot([eta, eta, xlim[0]], [ylim[0], 0, 0], "k")
ax[0,1].plot([mu, mu, xlim[0]], [ylim[0], 1, 1], 'k')

ax[0,1].plot(xlim, (xlim - eta)/(mu - eta), 'r')

# interpolations
r = np.array([-1.1, 0.3, 1.7])
N = (r-eta)/(mu-eta)

ax[0,1].fill_between([xlim[0], eta], [0,0], [ylim[0], ylim[0]],
                   color = colors[0], alpha = 0.8)

ax[0,1].fill_between([xlim[0], eta, eta, 1], [0,0,ylim[0], ylim[0]], [1,1,1,1],
                   color = colors[1], alpha = 0.5)
ax[0,1].fill_between([xlim[0], 1, 1, xlim[1]], [1,1,ylim[0], ylim[0]],
                   ylim[[1,1,1,1]],
                   color = colors[2], alpha = 0.5)



ax[0,1].set_ylim(ylim)
ax[0,1].set_xlim(xlim)

ax[0,1].set_xticks(np.append(r, [mu, eta]))
ax[0,1].set_xticklabels(["inter>intra", "inter<intra",
                         "sign(Inter)\n" r"$\neq$sign(Intra)",
                         "$\mu_i$","$\eta_i$"])

ax[0,1].get_xticklabels()[0].set_color(colors[0])
ax[0,1].get_xticklabels()[1].set_color(colors[1])
ax[0,1].get_xticklabels()[2].set_color(colors[2])


ax[0,1].set_ylim(ylim)
ax[0,1].set_xlim(xlim)

ax[0,1].set_yticks([0, 1])
ax[0,1].set_yticklabels([0, 1], fontsize = 16)

ax[0,0].set_title("A", loc = "left")
ax[0,1].set_title("B", loc = "left")

ax[0,1].set_xlabel("(hypothetical)\ninvasion growth rates", fontsize = 14)
#ax[1].set_ylabel("Niche differences\n$\mathcal{N}_i$")
ax[0,1].set_title("Interpretation of\nNiche differences")

###############################################################################
# fitness differences

#fig, ax = plt.subplots(1,2, figsize = (9,6))

zero = 0

F = (0-eta)/(mu - eta)
ylim = np.array([-0.1, 1.1])
xlim = np.array([-0.5, 1.1])

# interpolations
ax[1,0].plot([eta, eta, xlim[0]], [ylim[0], 0, 0], "k")
ax[1,0].plot([mu, mu, xlim[0]], [ylim[0], 1, 1], 'k')
ax[1,0].plot([zero, zero, xlim[0]], [ylim[0], F, F], 'k--')

ax[1,0].plot([eta, mu], [0, 1], 'r')

ax[1,0].set_ylim(ylim)
ax[1,0].set_xlim(xlim)

ax[1,0].set_xticks([eta, zero, mu])
ax[1,0].set_xticklabels(["$\eta_i$\n$\sim c_{ij}N_j^*$", "$0$\n$\sim N_i^*$", "$\mu_i$\n$\sim 0$"], fontsize = 16)

ax[1,0].set_yticks([0, F, 1])
ax[1,0].set_yticklabels([0, r"$\frac{0-\eta_i}{\mu_i-\eta_i}$", 1], fontsize = 16)

ax[1,0].set_xlabel("(hypothetical)\nequilibrium densities", fontsize = 14)
ax[1,0].set_ylabel("Fitness differences\n$\mathcal{F}_i$", fontsize = 14)
ax[1,0].set_title("Definition of\nFitness differences")

###############################################################################
# interpretation of fitness differences

colors = ["purple", "purple", "red"]

xlim = np.array([-1.8,2.4])
ylim = np.array([-1,2])

# interpolations
ax[1,1].plot([eta, eta, xlim[0]], [ylim[0], 0, 0], "k")
ax[1,1].plot([mu, mu, xlim[0]], [ylim[0], 1, 1], 'k')

ax[1,1].plot(xlim, (xlim - eta)/(mu - eta), 'r')

# interpolations
N = (r-eta)/(mu-eta)

ax[1,1].fill_between([xlim[0], eta], [0,0], [ylim[0], ylim[0]],
                   color = colors[0], alpha = 1)
ax[1,1].fill_between([xlim[0], eta], [0,0], [ylim[0], ylim[0]],
                   color = "k", alpha = 0.3)

ax[1,1].fill_between([xlim[0], eta, eta, 1], [0,0,ylim[0], ylim[0]], [1,1,1,1],
                   color = colors[1], alpha = 0.5)
ax[1,1].fill_between([xlim[0], 1, 1, xlim[1]], [1,1,ylim[0], ylim[0]],
                   ylim[[1,1,1,1]],
                   color = colors[2], alpha = 0.5)



ax[1,1].set_ylim(ylim)
ax[1,1].set_xlim(xlim)

ax[1,1].set_xticks(np.append(r, [mu, eta]))
ax[1,1].set_xticklabels(["higher\nequilibrium", "lower\nequilibrium", "negative\nequilibrium",  "$\mu_i$","$\eta_i$"])

ax[1,1].get_xticklabels()[0].set_color(colors[0])
ax[1,1].get_xticklabels()[1].set_color(colors[1])
ax[1,1].get_xticklabels()[2].set_color(colors[2])


ax[1,1].set_ylim(ylim)
ax[1,1].set_xlim(xlim)

ax[1,1].set_yticks([0, 1])
ax[1,1].set_yticklabels([0, 1], fontsize = 16)

ax[1,1].set_xlabel("(hypothetical)\nequilibrium densities", fontsize = 14)
#ax[1].set_ylabel("Niche differences\n$\mathcal{N}_i$")
ax[1,0].set_title("C", loc = "left")
ax[1,1].set_title("D", loc = "left")
ax[1,1].set_title("Interpretation of\nFitness differences")

fig.tight_layout()

fig.savefig("Figure_interpretation_NFD.pdf")