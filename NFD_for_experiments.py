"""
@author: J.W.Spaak
Numerically compute ND and FD for experimental data
"""

import numpy as np
import matplotlib.pyplot as plt
from warnings import warn

from scipy.interpolate import UnivariateSpline as uni_sp
from scipy.optimize import brentq

try:
    from numerical_NFD import NFD_model
except ImportError:
    # in case this code is used in a submodule, import from the submodule
    from nfd_definitions.numerical_NFD import NFD_model

class InputError(Exception):
    pass

def NFD_experiment(dens, time, r_i, N_star = None, na_action = "impute",
                 f0 = "spline", k = 3, s = None, log = True, visualize = True):
    """Compute the ND and FD for two-species experimental data
    
    Compute the niche difference (ND), niche overlapp (NO), 
    fitnes difference(FD) and conversion factors (c). The 3 experiments to
    conduct are described in:
    The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
    DOI:
    
    Parameters
    -----------
    N_star: None or ndarray (shape = 2)
        Monoculture equilibrium density for both species. If `None` N_star
        will be computed automatically.
    dens: ndarray (shape = (2, r, t))
        Densities of the species over time. dens[i,r,t] is the density of
        species i in replica r at time time[t]. The different starting
        densities can both be combined into different replicas.
        `np.nan` are allowed, but will be removed or imputed (see `na.action`).
        Alternatively (not recommended) `dens` can be a list-structure similar
        to (2,r,t) wherein each list can have separate length.
        If this is the case, `time` must also be a list, containing the
        datapoint for each timepoint.
    time: array like
        Timepoints at which measurments were taken in increasing order.
        If timepoints differ for different experiment `time` must contain all 
        timepoints. Alternatively, if `dens` is a list-structure time must have
        the same structure.
    na_action: "remove" or "impute" or "impute_geom" (default = "impute")
        If "remove" all datapoints with NA will be removed. This will cause the
        lost of the growth rate before and after the measured NA. Alternatively
        the missing value will be imputed using the datapoint before and after
        either using the mean ("impute") or the geometric mean ("impute_geom").        
    r_i: ndarray (shape = 2)
        Invasion growth rate of both species
    f0: "Linear", "spline" or ndarray, optional (shape = 2)
        Way to compute initial growth rate ``f0``. "Spline" will use the spline
        interpolation. "Linear" will assume constant growth between the first
        two data points. Alternatively ``f0`` can be passed directly.
    k : int or array of ints, optional, default = 3
        Degree of the smoothing spline.  Must be <= 5. If an array is passed
        the degree belongs to each experiment seperately in the following way
        [exp_1_sp_1, exp1_sp2, exp2_sp1, exp2_sp2]
    s : None, float or array of floats, optional
        Positive smoothing factor used to choose the number of knots. Number
        of knots will be increased until the smoothing condition is satisfied:
    
            sum((w[i] * (y[i]-spl(x[i])))**2, axis=0) <= s
    
        If None (default), ``s = len(w)`` which should be a good value if
        ``1/w[i]`` is an estimate of the standard deviation of ``y[i]``.
        If 0, spline will interpolate through all data points.
        If an array is passed the value belongs to each experiment seperately
        in the following way: [exp_1_sp_1, exp1_sp2, exp2_sp1, exp2_sp2]
    log: boolean, optinal, default = True
        If True fitting and visualizing is done one log transformed data
    visualize: boolean, optional, default = True
        If true, plot a graph of the growth rates and the percapita growth
        rates of both species. `fig`and `ax` of this figure are returned
        
    Returns
    -------
    pars : dict
        A dictionary with the following keys: 
            
    ``N_star`` : ndarray (shape = (n_spec, n_spec))
        N_star[i] equilibrium density with species `i`
        absent. N_star[i,i] is 0
    ``r_i`` : ndarray (shape = n_spec)
        invsaion growth rates of the species
    ``c`` : ndarray (shape = (n_spec, n_spec))
        The conversion factors from one species to the
        other. 
    ``ND`` : ndarray (shape = n_spec)
        Niche difference of the species to the other species
    ``NO`` : ndarray (shape = n_spec)
        Niche overlapp of the species (NO = 1-ND)
    ``FD`` : ndarray (shape = n_spec)
        Fitness difference
    ``f0``: ndarray (shape = n_spec)
        no-competition growth rate, f(0)
    fig: Matplotlib figure
        only returned if ``visualize`` is True
    ax: Matplotlib axes
        only returned if ``visualize`` is True
        
    Literature:
    The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
    DOI:
    """
    '''
    # check input experiment one must be increasing
    if (dens_exp1[:,:-1]>dens_exp1[:,1:]).any():
        warn("Densities in the first experiment are not strictly"
                         "increasing")
    if (dens_exp2[:,:-1]<dens_exp2[:,1:]).any():
        warn("Densities in the second experiment are not strictly"
                         "decreasing")'''            
    try:
        k[0]
    except TypeError:
        k = [k,k]
    try:
        s[0]
    except TypeError:
        s = [s,s]
    dens = np.array(dens)
    time = np.array(time)
    # per capita growth rate for both species in monoculture
    f, dict_N_t = per_capita_growth(dens, time, N_star,f0, k, s, log)


    if N_star is None:
        N_star = np.nanmean(dens, axis = (1,2))
    # compute the ND, FD etc. parameters
    pars = {"N_star": np.array([[0,N_star[1]],[N_star[0],0]]),
            "r_i": r_i}    
    pars = NFD_model(f, pars = pars, experimental = True)
    if visualize: # visualize results if necessary
        fig, ax = visualize_fun(f,times, exps, N_star, pars, dict_N_t, log)
        return pars, dict_N_t, fig, ax
    
    return pars, dict_N_t
    
def per_capita_growth(dens, time, N_star, f0, k, s, log):
    """interpolate the per capita growth rate of the species
    
    times:
        Timepoints of measurments
    exps:
        Densities of the species at timepoints times
    N_star: float
        equilibrium density of the species
    f0: float
        monoculture growth rate
    k: int
        Order of interpolation spline
    s: float or None
        Smoothing factor
    log: boolean, optinal, default = True
        If True fitting and visualizing is done one log transformed data
        
    Returns
    -------
    f: callable
        Per capita growth rate of species, fullfilling the differential
        equation dN/dt=N*f(N)
        Values below min(exps) are assumed to be f0, Values above max(exps)
        are assumed to be f(max(exps))
    """
    # percapita growth rates for each of the experiments separately
    dict_subf = {} # contains f for each sub experiment and species
    dict_N_t = {} # contains fitted growth curves for each sub experiment
    for i in range(2):
        f_N, N_t = dens_to_per_capita(dens[i], time, k[i], s[i], log, f0[i])
        dict_subf["f_spec{}".format(i)] = f_N
        dict_N_t["spec{}".format(i)] = N_t
    
    # monoculture growth rate, assumes that growth rate is constant at
    # the beginning
    if f0 == "spline":
        f0 = [dict_subf["f_spec"+str(i)](np.nanmin(dens[i])) for i in [0,1]]
        print("f0", f0)
    elif f0 == "linear":
        f0 = np.log(exps[1][:,1]/exps[1][:,0])/(times[1][1]-times[1][0])
    # other wise f0 is assumed to be an array of shape (2,)
            
    
    # per capita growth rate for each species, using different cases
    def f_spec(N,i):
        # `N` species density, `i` speces index
        if N<np.nanmin(dens[i]): # below minimum, use f0
            return f0[i]
        elif N>np.nanmax(dens[i]): # above maximum, use maximal entry
            return dict_subf["f_spec"+str(i)](np.nanmax(dens[i]))
        else: # above maximum
            return dict_subf["f_spec"+str(i)](N)
    
    def f(N):
        """ per capita growth rate of species
        
        can only be used for monoculture, i.e. f(N,0) or f(0,N).
        growth rate of non-focal species is set to np.nan"""
        if np.all(N==0):
            return np.array([f_spec(0,0), f_spec(0,1)])
        else:
            spec_foc = np.argmax(N)
        ret = np.full(2, np.nan)
        ret[spec_foc] = f_spec(max(N),spec_foc)
        return ret
    
    return f, dict_N_t
            
def dens_to_per_capita(dens, time, k, s, log, f0 = None):
    # convert densities over time to per capita growth rate
    time_diff = time[1:]-time[:-1]
    per_cap_growth = np.log(dens[...,1:]/dens[...,:-1])/time_diff

    # remove nan's
    finite = np.isfinite(dens[...,:-1]) & np.isfinite(per_cap_growth)
    dens_finite = dens[...,:-1][finite]
    per_cap_growth = per_cap_growth[finite]
    # sort for increasing dens, needed for spline
    ind = np.argsort(dens_finite)
    dens_finite = dens_finite[ind]
    per_cap_growth = per_cap_growth[ind]
    if not (f0 is None):
        per_cap_growth[0] = f0
        w = np.ones(len(dens_finite))
        w[0] = 1e10 # to force spline through first point
    else:
        w = np.ones(len(dens_finite))
    
    dNdt = uni_sp(np.log(dens_finite), per_cap_growth, k = k, s = s, w = w)
    def per_capita(N):
        # compute the percapita growth rate of the species        
        return dNdt(np.log(N))
    N_t = None
    return per_capita, N_t

def visualize_fun(f,times, exps, N_star, pars, dict_N_t, log):
    # visualize the fitted population density and the per capita growthrate
    fig, ax = plt.subplots(2,2, figsize = (11,11),
                           sharey = "row", sharex ="row")
    
    # plot the densities over time
    # plot real data of exp1 for both species
    ax[0,0].scatter(times[1], exps[1][0], color = "black")
    ax[0,1].scatter(times[1], exps[1][1], color = "black")
    
    # plot real data of exp2 for both species
    ax[0,0].scatter(times[2], exps[2][0],facecolor = "none", color = "black")
    ax[0,1].scatter(times[2], exps[2][1],facecolor = "none", color = "black")
    
    # time lines
    time_1 = np.linspace(*times[1][[0,-1]], 100)
    time_2 = np.linspace(*times[2][[0,-1]], 100)
    
    # plot fitted data of exp1 for first
    ax[0,0].plot(time_1, dict_N_t["exp1_spec0"](time_1),label = "fit, exp1")
    ax[0,0].plot(time_2, dict_N_t["exp2_spec0"](time_2),label = "fit, exp2")
    
    # plot fitted data of exp1 for second
    ax[0,1].plot(time_1, dict_N_t["exp1_spec1"](time_1),label = "fit, exp1")
    ax[0,1].plot(time_2, dict_N_t["exp2_spec1"](time_2),label = "fit, exp2")
    
    ax[0,0].axhline(N_star[0], linestyle = "dotted", color = "black")
    ax[0,1].axhline(N_star[1], linestyle = "dotted", color = "black")
    
    ax[0,0].legend()
    ax[0,1].legend()
    
    # add axis labeling   
    ax[0,0].set_title("Species 1, densities")
    ax[0,1].set_title("Species 2, densities")
    
    ax[0,0].set_ylabel(r"Densities $N_i(t)$")
    ax[0,0].set_xlabel(r"Time $t$")
    ax[0,1].set_xlabel(r"Time $t$")
    
    # log transform y axis
    if log:
        ax[0,0].semilogy()
        ax[0,1].semilogy()
    
    # plot the fitted per capita growth rate
    if log:
        N_1 = np.exp(np.linspace(np.log(exps[1][0,0]), np.log(exps[2][0,0])
                ,100))
        N_2 = np.exp(np.linspace(np.log(exps[1][1,0]), np.log(exps[2][1,0])
                ,100))
        ax[1,0].semilogx()
        ax[1,1].semilogx()
    else:
        N_1 = np.linspace(exps[1][0,0], exps[2][0,0],100)
        N_2 = np.linspace(exps[1][1,0], exps[2][1,0],100)
    
    ax[1,0].plot(N_1,[f([N,0])[0] for N in N_1])
    ax[1,1].plot(N_2,[f([0,N])[1] for N in N_2])
    
    # add equilirium and 0 axis line
    ax[1,0].axhline(0,linestyle = "dotted", color = "black")
    ax[1,0].axvline(N_star[0],linestyle = "dotted", color = "black")
    ax[1,0].text(N_star[0], ax[1,0].get_ylim()[1]*0.9,r"equi. $N_1^*$"
          , rotation = 90, fontsize = 14, ha = "right")
    
    ax[1,1].axhline(0,linestyle = "dotted", color = "black")
    ax[1,1].axvline(N_star[1],linestyle = "dotted", color = "black")
    ax[1,1].text(N_star[1], ax[1,1].get_ylim()[1]*0.9,r"equi. $N_2^*$"
          , rotation = 90, fontsize = 14, ha = "right")
    
    # add points assuming constant growth between two measurments
    growth_1 = np.log(exps[1][:,1:]/exps[1][:,:-1])
    growth_2 = np.log(exps[2][:,1:]/exps[2][:,:-1])
    per_cap_1 = growth_1/(times[1][1:] - times[1][:-1])
    per_cap_2 = growth_2/(times[2][1:] - times[2][:-1])
    
    if log: # take geometric average between two densities
        av_dens_1 = np.sqrt(exps[1][:,1:] * exps[1][:,:-1])
        av_dens_2 = np.sqrt(exps[2][:,1:] * exps[2][:,:-1])
    else:
        av_dens_1 = (exps[1][:,1:] + exps[1][:,:-1])/2
        av_dens_2 = (exps[2][:,1:] + exps[2][:,:-1])/2
    
    # plot the measured per capita growth rates
    ax[1,0].plot(av_dens_1[0], per_cap_1[0], 'o', label = "measured")
    ax[1,0].plot(av_dens_2[0], per_cap_2[0], 'o')
    ax[1,1].plot(av_dens_1[1], per_cap_1[1], 'o')
    ax[1,1].plot(av_dens_2[1], per_cap_2[1], 'o')
    
    # add point where ND equality was computed
    x_dist = ax[1,0].get_xlim()[1]-ax[1,0].get_xlim()[0]
    c_N_star = pars["c"][[0,1],[1,0]]*N_star[[1,0]]
    ax[1,0].plot(c_N_star[0], f([c_N_star[0],0])[0], 'o')
    
    ax[1,0].text(c_N_star[0]+0.03*x_dist, f([c_N_star[0],0])[0],
                  r"$c_2\cdot N_2^*$", fontsize =14)
    ax[1,1].plot(c_N_star[1], f([0,c_N_star[1]])[1], 'o')
    ax[1,1].text(c_N_star[1]+0.03*x_dist, f([0,c_N_star[1]])[1],
                  r"$c_1\cdot N_1^*$", fontsize = 14)
    
    # axis labeling
    ax[1,0].set_title("Species 1, per capita growth")
    ax[1,1].set_title("Species 2, per capita growth")
    
    ax[1,0].set_ylabel(r"Percapita growth rate $f_i(N_i,0)$")
    ax[1,0].set_xlabel(r" Density $N_1$")
    ax[1,1].set_xlabel(r" Density $N_2$")
    return fig, ax