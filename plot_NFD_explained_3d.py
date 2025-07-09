import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize = (9,4.5))

###############################################################################
# coexistence interpretation of niche and fitness differences
ax = fig.add_axes([.1,.2, 0.35, .7])
zero = 0
mu = 1
eta = -0.4
r = 0.6

N = (r-eta)/(mu - eta)
F = (0-eta)/(mu-eta)

ylim = np.array([-0.1, 1.1])
xlim = np.array([-0.5, 1.1])

# interpolations
ax.plot([eta, eta, xlim[0]], [ylim[0], 0, 0], "k")
ax.plot([mu, mu, xlim[0]], [ylim[0], 1, 1], 'k')
ax.plot([zero, zero, xlim[0]], [ylim[0], F, F], 'k--')
ax.plot([r, r, xlim[0]], [ylim[0], N, N], 'k--')

ax.plot([eta, mu], [0, 1], 'k')

ax.set_ylim(ylim)
ax.set_xlim(xlim)

ax.set_xticks([eta, zero,r,  mu])
ax.set_xticklabels(["$\eta_i$", "$0$", "$r_i$", "$\mu_i$"], fontsize = 16,
                   color = "b")
ax.set_xlabel("Growth rates", fontsize = 16)

ax.set_yticks([0, F, N, 1])
ax.set_yticklabels([0, r"$\mathcal{F}_i=$" + "\n" + r"$\frac{0-\eta_i}{\mu_i-\eta_i}$",
                    r"$\mathcal{N}_i=$" + "\n" + r"$\frac{r_i-\eta_i}{\mu_i-\eta_i}$",1],
                   fontsize = 16, color = "g")

ax.set_ylabel("Niche and Fitness differences", fontsize = 14)
ax.set_title("Interpretation of coexistence")
###############################################################################
# 3d plot
ax = fig.add_axes([.41,0, 0.6, 1], projection = "3d")
ax.set_axis_off()

itera = 200
zeros = np.zeros(itera)
x = np.linspace(0,1.9,itera)
ax.plot(x,zeros, zeros, 'k')
ax.plot(zeros,x,zeros, 'k')
ax.plot(zeros, zeros, np.linspace(0, 1.2, itera), 'k')
ax.plot(zeros, zeros, np.linspace(0,-0.65, itera), 'k--')

ax.set_xlim(x[[0,-1]])
ax.set_ylim(x[[0,-1]])
ax.set_zlim(0.2, 1.3)

ax.view_init(ax.elev, ax.azim+90)
fs = 14
col_growth = "b"
col_dens = "r"
col_NFD = "green"
col_conn = "orange"

# growth rates
mu_i = 1

# relevant densities
N_j_star = 1
c_ijNj = 1
N_i_star = 0.62
N_j_dagger = 1.8 # where species 1 has no growth rate

# growth rates 
r_i = mu_i*(N_j_star - N_j_dagger)/(0-N_j_dagger)
eta_i = mu_i*(c_ijNj - N_i_star)/(0-N_i_star)

# plot growth rates
ax.plot(0,0,mu_i, '.', color = col_dens)
ax.text(0.05,-0.04,mu_i, r'$\mu_i$', va = "bottom", ha = "right", fontsize = fs,
        color = col_growth)

ax.plot(N_i_star, 0, 0, '.', color = col_dens)
ax.text(N_i_star, 0,0, '$N_i^*$', color = col_dens,
        va = "bottom", ha = "right", fontsize = fs)

ax.plot(c_ijNj, 0,0, '.', color = col_dens)
ax.text(c_ijNj, 0,0, "$c_{ij}N_j^*$", color = col_dens,
        va = "bottom", ha = "right", fontsize = fs)

ax.plot(0, N_j_star, 0, '.', color = col_dens)
ax.text(0, N_j_star, 0+0.09, " $N_j^*$", color = col_dens,
        va = "center", ha = "left", fontsize = fs)

ax.plot(0, N_j_dagger, 0, '.', color = col_dens)
ax.text(0, N_j_dagger, 0+0.09, " $N_j^\dagger$", color = col_dens,
        ha = "left", va = "center", fontsize = fs)

ax.plot([N_i_star, 0], [0, N_j_dagger], [0,0], '--', color = col_growth,
        label = "$f_i(N_i, N_j)=0$")

# plot growth rates plane
normal_vec = np.cross([-N_i_star,0, 1], [0, -N_j_dagger, 1])
xx, yy = np.meshgrid(x, x)
offset = normal_vec[2]*mu_i
zz = (offset - normal_vec[0]*xx - normal_vec[1]*yy)/normal_vec[2]
zz[zz<-1.3] = np.nan
surf = ax.plot_surface(xx, yy, zz, alpha = 0.2, color = col_growth, label = "$f_i(N_i, N_j)$")
surf._edgecolors2d = surf._edgecolor3d
surf._facecolors2d = surf._facecolor3d
ax.legend()

# plot growth rates
ax.plot(2*[c_ijNj], [0,0], [0, eta_i], ':', color = col_conn)
ax.plot(c_ijNj, 0, eta_i, '.', color = col_growth)
ax.text(c_ijNj, 0, eta_i, '$\eta_i$', color = col_growth,
        va = "top", ha = "center", fontsize = fs)


# plot growth rates
ax.plot(2*[0], 2*[N_j_star], [0, r_i], ':', color = col_conn)
ax.plot(0, N_j_star, r_i, '.', color = col_growth)
ax.text(0, N_j_star, r_i, '$r_i$', color = col_growth,
        va = "bottom", ha = "left", fontsize = fs)

# add niche and fitness differences information
ax.plot([c_ijNj, 0], 2*[0], 2*[eta_i], ':', color = col_conn)
ax.plot(0,0,eta_i, '.', color = col_NFD)
ax.text(0,0, eta_i, "$0$", color = col_NFD, fontsize = fs)

ax.plot(2*[0], [N_j_star, 0], 2*[r_i], ':', color = col_conn)
ax.plot(0,0, r_i, '.', color = col_NFD)
ax.text(0,0, r_i, "$\mathcal{N}_i$", color = col_NFD, fontsize = fs)

ax.plot([N_i_star, 0], 2*[0], 2*[0], ":", color = col_conn)
ax.text(0,0,0, '$\mathcal{F}_i$', color = col_NFD, fontsize = fs)
ax.plot(0,0,0, '.', color = col_NFD)

ax.plot(0,0,mu_i, '.', color = col_NFD)
ax.text(-0.04,0.04,mu_i, '$1$', color = col_NFD, ha = "left", va = "bottom",
        fontsize = fs)
ax.text(0,0,mu_i+0.09, '$\mapsto$', color = "k", ha = "center", va = "bottom",
        fontsize = fs)

# compare niche and fitness differences
N_i_hash = N_i_star/c_ijNj * (N_j_star  - N_j_dagger)/(0-N_j_dagger)
ax.plot([0, N_i_hash], 2*[N_j_star], 2*[0], ':', color = col_conn)
ax.plot(N_i_hash, N_j_star, 0, '.', color = col_dens)
ax.text(N_i_hash, N_j_star, 0, "$\mathcal{N}_i-\mathcal{F}_i$", color = col_NFD,
        ha = "center", va = "top", fontsize = fs)
ax.plot(2*[N_i_hash], [0,N_j_star], 2*[0], ":", color = col_conn)
ax.plot(N_i_hash, 0, 0, '.', color = col_dens)
ax.text(N_i_hash, 0, 0, "$N_i^\#$", color = col_dens,
        ha = "right", fontsize = fs, va = "bottom")

fig.savefig("Figure_explain_NFD.pdf")