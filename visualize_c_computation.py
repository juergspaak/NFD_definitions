import numpy as np
import matplotlib.pyplot as plt


import numerical_NFD_permanence as NF


from matplotlib import cm

def plot_c_computations(pars, c_range = None, sharey = True,
                        names = None, colors = None):
    n_spec = len(pars["r_i"])
    
    if c_range is None:
        c_save = pars["c"].copy()
        c_save[c_save == 0] = 1e-5
        c_save[np.isinf(c_save)] = 1e5
        c_range = np.geomspace(np.nanmin(c_save)/5, np.nanmax(c_save)*5,
                               101)
    if isinstance(c_range, int):
        c_save = pars["c"].copy()
        c_save[c_save == 0] = 1e-5
        c_save[np.isinf(c_save)] = 1e5
        c_range = np.geomspace(np.nanmin(c_save)/5, np.nanmax(c_save)*5,
                               c_range)
    if names is None:
        names = np.arange(n_spec)
    fig, ax = plt.subplots(n_spec, n_spec, figsize = (3*n_spec, 3*n_spec),
                           sharex = True, sharey = sharey)
    ax[0,0].set_xlim(c_range[[0,-1]])
    # extend c_range to include all densities in the range c and 1/c
    c_range = np.concatenate([1/c_range[1/c_range>np.amax(c_range)],
                              c_range[::-1],
                             1/c_range[1/c_range<np.amin(c_range)]])[::-1]
    if colors is None:
        colors = cm.viridis(np.linspace(0,1, n_spec))
    for i in range(n_spec):
        for j in range(n_spec):
            
            NO_ij = NF.NO_fun(pars, [i,j])
            
            # compute the niche overlapp with species i as a function of c
            NO_ij = [np.abs(NO_ij(c)) for c in c_range]
            
            
            
            if i!= j:
                ax[i,j].axvline(pars["c"][i,j], color = "k", label = "c_ij")
                ax[i,j].semilogx(c_range, NO_ij,
                                 color = colors[i])
                ax[j,i].semilogx(1/c_range, NO_ij, color = colors[i])
            
            
    for i in range(n_spec):
        ax[i,i].text(0.5, 0.5, names[i],
                     color = colors[i], fontsize = 18,
                     va = "center", ha = "center",
                     transform=ax[i,i].transAxes)
    
        
    return fig, ax