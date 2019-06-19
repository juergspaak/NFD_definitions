"""
@author: J.W.Spaak
Numerically compute ND and FD for experimental data
"""

import numpy as np
import matplotlib.pyplot as plt
from warnings import warn

from scipy.interpolate import UnivariateSpline as uni_sp
from scipy.optimize import brentq, fsolve
from scipy.integrate import odeint

try:
    from numerical_NFD import NFD_model
except ImportError:
    # in case this code is used in a submodule, import from the submodule
    from nfd_definitions.numerical_NFD import NFD_model

class InputError(Exception):
    pass

def NFD_experiment(dens, time, r_i, N_star = None, na_action = "remove",
                 f0 = "spline", k = 3, s = "fac=1", log = True,
                 id_exp_1 = None, visualize = True, extrapolate = "True",
                 growth_data = np.zeros((2,2,0))):
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
    na_action: "remove" or "impute" or "impute_geom" (default = "remove")
        If "remove" all datapoints with NA will be removed. This will cause the
        lost of the growth rate before and after the measured NA. Alternatively
        the missing value will be imputed using a running mean ("impute") or
        a running geometric mean ("impute_geom") on 3 datapoints.
        Imputing NA is assuming equidistant timepoints, if this is not the case
        imputing NA is not recommended.
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
        monoculture growth rate, f(0)
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
    input_par, log, id_exp_1 = __input_check_exp__(k , s, f0, N_star,
                    id_exp_1, log, extrapolate, dens, time, growth_data)         
    
    dens = np.array(dens)
    time = np.array(time)
        
    # impute na if necessary
    if na_action == "impute": # compute mean
        imp_av = lambda X, axis: log["exp"](np.nanmean(log["log"](X),
                                    axis = axis))
        na_imputed = dens.copy()
        na_imputed[...,1:-1] = imp_av([na_imputed[...,:-2],
                na_imputed[...,1:-1], na_imputed[...,2:]], axis = 0)
        dens[np.isnan(dens)] = na_imputed[np.isnan(dens)] # replace imputation
    
    # per capita growth rate for both species in monoculture
    f, f_spec = per_capita_growth(dens, time, input_par, log)

    if N_star is None:
        N_star = log["nanmean"](dens[...,-1], axis = (1)) # starting guess
        # compute N_star
        for i in range(2):
            N_star[i] = fsolve(f_spec, N_star[i], args = (i,))
    c = np.ones((2,2))
    c[[1,0],[0,1]] = [N_star[0]/N_star[1], N_star[1]/N_star[0]]    
    # compute the ND, FD etc. parameters
    pars = {"N_star": np.array([[0,N_star[1]],[N_star[0],0]]),
            "r_i": r_i, "f": f, "c": c}
    pars = NFD_model(f, pars = pars, experimental = True)

    N_t_fun, N_t_data = dens_over_time(f_spec, time, dens, id_exp_1)
    if visualize: # visualize results if necessary
        fig, ax = visualize_fun(time, dens, pars, N_t_fun, id_exp_1, log, f_spec)
        return pars, N_t_fun, fig, ax
    
    return pars

def __input_check_exp__(k , s, f0, N_star, id_exp_1, log, extrapolate, 
                        dens, time, growth_data):
    # convert input parameters to right format if necessary
    input_par = {"k":k, "s":s}
    for key in input_par.keys():
        if type(input_par[key]) != list:
            input_par[key] = 2*[input_par[key]]
            
    input_par["extrapolate"] = extrapolate   
 
    if id_exp_1 is None:
        id_exp_1 = np.full(dens.shape[1], False)
        id_exp_1[:len(id_exp_1)//2] = True
        
    if log: # fitting, imputing, plotting etc is all done in logspace
        geo_nanmean = lambda x, axis = None: np.exp(np.nanmean(np.log(x)
                , axis = axis))
        geo_mean = lambda x, axis = None: np.exp(np.mean(np.log(x)
                , axis = axis))
        log = {"log": np.log, "exp": np.exp, "log_space": True,
               "nanmean": geo_nanmean, "mean": geo_mean}
    else: # log and exp functions are identity functions
        log = {"log": lambda x: x , "exp": lambda x: x, "log_space": False,
               "nanmean": np.nanmean, "mean": np.mean}
        
    # compute data intern variance
    time_diff = time[1:]-time[:-1]
    per_cap_growth = np.log(dens[...,1:]/dens[...,:-1])/time_diff
    input_par["var"] = np.nanmean(np.nanvar(per_cap_growth, axis = 1),
             axis = -1)
    
    # add predefined N_star value to growth_data if needed
    if not (N_star is None):
        # per cap grwth is 0 at equilibrium
        growth_append = np.array([N_star, [0,0]]).T.reshape(2,2,1)
        growth_data = np.append(growth_data, growth_append, axis = -1)
            
    
    
    if (type(f0) == str) and f0[:5] == "perc=":
        f0 = np.nanpercentile(per_cap_growth, float(f0[5:]), axis = (1,2))
    elif f0 == "linear":
        f0 = np.nanmean(per_cap_growth[:,id_exp_1,0], axis = 1)
    print(f0)    
    if not (f0 is None):
        scal = 5.0**(-np.arange(4).reshape(-1,1))
        dens_ap = scal*np.nanmin(dens, axis = (1,2))
        gr_ap = f0*np.ones((4,2))
        growth_append = [[dens_ap[:,0], gr_ap[:,0]],[dens_ap[:,1], gr_ap[:,1]]]
        growth_data = np.append(growth_data, growth_append, axis = -1)        
        
    input_par["growth_data"] = growth_data
    return input_par, log, id_exp_1
    
    
def per_capita_growth(dens, time, input_par, log):
    """interpolate the per capita growth rate of the species
    
    times:
        Timepoints of measurments
    exps:
        Densities of the species at timepoints times
    N_star: float
        equilibrium density of the species

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
    dict_f = {} # contains per capita growth rate for each specie in monocul.
    for i in range(2):
        f_N = dens_to_per_capita(dens[i], time, input_par, log, i)
        dict_f["spec{}".format(i)] = f_N     
    
    # per capita growth rate for each species, using different cases
    mins = np.nanmin(dens, axis = (1,2))
    maxs = np.nanmax(dens, axis = (1,2))
    
    def f_spec(N,i):
        # `N` species density, `i` speces index
        # below minimum, use f0
        if (N < mins[i]): 
            return dict_f["spec{}".format(i)](mins[i])
        # above maximum, use maximal entry
        elif (N > maxs[i]) and not input_par["extrapolate"]: 
            return dict_f["spec{}".format(i)](maxs[i])
        else: # above maximum
            return dict_f["spec{}".format(i)](N)
    
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
    
    return f, f_spec

def dens_to_per_capita(dens, time, input_par, log, i):
    # convert densities over time to per capita growth rate
    time_diff = time[1:]-time[:-1]
    per_cap_growth = np.log(dens[...,1:]/dens[...,:-1])/time_diff

    # use middle point for growth rate, similar in plotting function
    dens_mid = log["nanmean"]([dens[...,1:], dens[...,:-1]], axis = 0)
    
    # remove nan's
    finite = np.isfinite(dens_mid) & np.isfinite(per_cap_growth)
    dens_finite = dens_mid[finite]
    per_cap_growth = per_cap_growth[finite]
    w = np.ones(len(dens_finite)) # weighting for univariate spline
    
    # add predefined growth data
    print("here", input_par["growth_data"][i,0])
    print(input_par["growth_data"])
    dens_finite = np.append(dens_finite, input_par["growth_data"][i,0])
    per_cap_growth = np.append(per_cap_growth, input_par["growth_data"][i,1])
    w = np.append(w, np.full(len(dens_finite)-len(w),input_par["var"][i]*1e10))
    # sort for increasing dens, needed for spline
    ind = np.argsort(dens_finite)
    dens_finite = dens_finite[ind]
    per_cap_growth = per_cap_growth[ind]
    w = w[ind]
    
    # remove double values (can occur from imputing na)
    unique = dens_finite[1:] > dens_finite[:-1]
    unique = np.append([True], unique)
    dens_finite = dens_finite[unique]
    per_cap_growth = per_cap_growth[unique]
    w = w[unique]
    
    # compute smoothingfactor s
    if type(input_par["s"][i]) == float:
        s = input_par["s"][i]
    elif ((type(input_par["s"][i]) == str) and 
                (input_par["s"][i][:4] == "fac=")):
        s = float(input_par["s"][i][4:])*input_par["var"][i]*len(w)
    else:
        s = None
    dNdt = uni_sp(log["log"](dens_finite), per_cap_growth,
                  k = input_par["k"][i], s = s, w = w)
    
    def per_capita(N):
        # compute the percapita growth rate of the species        
        return dNdt(log["log"](N))
    
    return per_capita

def dens_over_time(f_spec, time, dens, id_exp_1):
    # compute the densities over time for both species
    N_t_fun = {} # contains a function to compute densities over time
    N_t_data = {} # contains the densities over time
    for i in range(2):
        for start in ["low", "high"]:
            key = "spec{}_{}".format(i, start)
            if start == "low":
                N_start = np.nanmean(dens[i, id_exp_1, 0])
            else:
                N_start = np.nanmean(dens[i, ~id_exp_1, 0])
            N_t_fun[key] = lambda time, i = i, N_start = N_start: odeint(
                    lambda N,t: N*f_spec(N,i), N_start, np.append(0,time))[1:]
            N_t_data[key] = N_t_fun[key](time)
    return N_t_fun, N_t_data
    
def visualize_fun(time, dens, pars, N_t_fun, id_exp_1, log, f_spec):
    # visualize the fitted population density and the per capita growthrate
    fig, ax = plt.subplots(2,2, figsize = (12,12),
                           sharey = "row", sharex ="row")
    # plot the measured densities and the fitted densities
    col = ["black", "grey"]
    for i in range(2):
        j = [1,0][i] # the index of the other species
    
        # plot the densities over time
        # plot real data of exp1 for both species
        ax[0,i].plot(time, dens[i, id_exp_1].T, '^', color = col[i])
        
        # plot real data of exp2 for both species
        ax[0,i].plot(time, dens[i, ~id_exp_1].T, 'o', color = col[i])
        
        # time lines
        time_fine = np.linspace(*time[[0,-1]], 100)
        # plot fitted data of exp1 and exp2
        ax[0,i].plot(time_fine, N_t_fun["spec{}_low".format(i)](time_fine),
          '-', color = "black", label = "fit, exp_low")
        ax[0,i].plot(time_fine, N_t_fun["spec{}_high".format(i)](time_fine),
          '--', color = "black", label = "fit, exp_high")
        
        ax[0,i].axhline(pars["N_star"][j,i], color = "green")
        
        # add axis labeling   
        ax[0,i].set_title("A; Species {}, densities".format(i+1))
        
        ax[0,0].set_ylabel(r"Densities $N_i(t)$")
        ax[0,i].set_xlabel(r"Time $t$")
        ax[0,i].legend()
        

    # plot the per capita growth rate
    time_diff = time[1:]-time[:-1]
    per_cap_growth = np.log(dens[...,1:]/dens[...,:-1])/time_diff
    for i in range(2):
        j = [1,0][i] 
        ax[1,i].plot(dens[i, :, :-1], per_cap_growth[i], 'o', color = col[i])
        dens_range = log["exp"](np.linspace(
                *np.nanpercentile(log["log"](dens[i]),[0,100]), 100))
        ax[1,i].plot(dens_range, [f_spec(de, i) for de in dens_range],
          color = "blue", label = r"$f_{}(N,0)$".format(i+1))
        
        # add equilibrium densities
        ax[1,i].axhline(0, linestyle = '--', color = "grey")
        ax[1,i].axvline(pars["N_star"][j,i], linestyle = "-", color = "green",
          label = r"$N^*_{}$".format(i+1))
        ax[1,i].axhline(pars["r_i"][i], linestyle = ':', color = "orange",
          label = r"$r_{}$ (invasion)".format(i+1))
         
    if log["log"]:
        # change to log scale
        ax[0,0].semilogy()
        ax[1,1].semilogx()
    
    # add converted densities of other species
    for i in range(2):
        j = [1,0][i] 
        N_star_j = (pars["N_star"]*pars["c"])[i,j]
        exp_id = "low" if N_star_j < pars["N_star"][j,i] else "high"
        
        ax[1,i].plot(N_star_j, f_spec(N_star_j, i), 'ro',
          label = r"$c_{}N_{}^*$".format(j+1,j+1))
        
        if N_star_j > np.nanmax(dens[i]):
            extend = log["exp"](np.linspace(log["log"](np.nanmax(dens[i])),
                        log["log"](N_star_j)))
            ax[1,i].plot(extend, [f_spec(N, i) for N in extend],
              "--", color = "blue")
            continue
        elif N_star_j < np.nanmin(dens[i]):
            extend = log["exp"](np.linspace(log["log"](np.nanmin(dens[i])),
                        log["log"](N_star_j)))
            ax[1,i].plot(extend, [f_spec(N, i) for N in extend],
              "--", color = "blue")
        
        # find time where density equates N_star_j
        N_t = N_t_fun["spec{}_{}".format(i, exp_id)]
        try:
            t_N_star_j = brentq(lambda t: (N_t(t) - N_star_j), 0, time[-1])
            ax[0,i].plot(t_N_star_j, N_t(t_N_star_j), 'ro',
              label = r"$c_{}N_{}^*$".format(j+1,j+1))
            
        except ValueError: # growth did not reach N_star_j in time
            ax[0,i].plot(time[-1], N_star_j, 'ro',
              label = r"$c_{}N_{}^*$".format(j+1,j+1))
    ax[1,1].legend(fontsize = 14)
    ax[1,0].legend(fontsize = 14)

    return fig, ax
fmin =  np.log(0.9)/3.5
pars, N_t_fun, fig, ax = NFD_experiment(dens, mono_days, r_i, visualize = True, 
                                   s = "fac=2", f0 = "perc=100", N_star = [5e5, 5e5],
                                   log = True)