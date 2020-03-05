"""
@author: J.W.Spaak
Numerically compute ND and FD for experimental data
"""
import numpy as np

from scipy.interpolate import UnivariateSpline as uni_sp
from scipy.optimize import brentq, fsolve
from scipy.integrate import odeint

try:
    from numerical_NFD import NFD_model, InputError
except ImportError:
    # in case this code is used in a submodule, import from the submodule
    from nfd_definitions.numerical_NFD import NFD_model, InputError
 
def NFD_experiment(dens, time, r_i, N_star = "average", na_action = "remove",
                 f0 = "spline", k = 3, s = "fac=1", log = True,
                 id_exp_1 = None, visualize = True, extrapolate = "True",
                 growth_data = np.zeros((2,2,0)), from_R = False):
    """Compute the ND and FD for two-species experimental data
    
    Compute the niche difference (ND), niche overlapp (NO), 
    fitnes difference(FD) and conversion factors (c). The 3 experiments to
    conduct are described in:
    The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
    DOI:
    
    Parameters
    -----------
    dens: ndarray (shape = (2, r, t))
        Densities of the species over time. dens[i,r,t] is the density of
        species `i` in replica r at time time[t]. The different starting
        densities should both be combined into different replicas.
        `np.nan` are allowed, but will be removed or imputed (see `na.action`).
    time: array like
        Timepoints at which measurments were taken in increasing order.
        If timepoints differ for different experiment `time` must contain all 
        timepoints.
    r_i: ndarray (shape = 2)
        Invasion growth rate of both species
    N_star: "spline", "average" or ndarray (shape = 2)
        Monoculture equilibrium density for both species. If "average" N_star
        will be set to the average at the last time point. If "splin" N_star
        will be computed using spline interpolation.
        Default is "average".
    na_action: "remove" or "impute"(default = "remove")
        If "remove" all datapoints with NA will be removed. This will cause the
        loss of the growth rate before and after the measured NA.
        Alternatively ("impute"), the missing value will be imputed using a
        running mean on 3 datapoints. Imputing NA is assuming equidistant
        timepoints, if this is not the case imputing NA is not recommended.
    f0: "linear", "perc=x", "spline" or ndarray (shape = 2), default = "linear"
        Way to compute initial growth rate ``f0``. "Spline" will use the spline
        interpolation. "Linear" will assume constant growth between the first
        two data points. "perc=x" will set f0 to the x percentile of the
        measured per-capita growth rates of said species.
        Alternatively ``f0`` can be passed directly.
    k : int or array of ints, optional, default = 3
        Degree of the smoothing spline.  Must be <= 5. If an array is passed
        the degree belongs to each species seperately in the following way
        [k_sp_1, k_sp2]
    s : None, float or "fac=x", where x is a float, or a list of those
        Positive smoothing factor used to choose the number of knots. Number
        of knots will be increased until the smoothing condition is satisfied:
    
            sum((w[i] * (y[i]-spl(x[i])))**2, axis=0) <= s
        
        Where y is the per_capita growth rate and x is the density.
        If 0, spline will interpolate through all data points.
        If None (default), ``s = len(y)`` which should be a good value if
        ``1`` is an estimate of the standard deviation of ``y[i]``.
        If "fac=x" then s=x*var, where var is the standard deviation of the
        per capita growth rates. Default is "fac=1"        
        If a list is passed the value belongs to each species seperately.
    log: boolean, optinal, default = True
        If True fitting and visualizing is done one log transformed data
    extrapolate: boolean, default = True
        If True, data above the measured densities will be extrapolated,
        otherwise they are set to the boundary values
    id_exp_1: boolean array or None
        Identity of replicas with low starting densities. If None, 
        the first half ot the replicas is assumed to be with low starting
        density.
    visualize: boolean, optional, default = True
        If true, plot a graph of the growth rates and the percapita growth
        rates of both species. `fig`and `ax` of this figure are returned
    growth_data: array of shape (2,2,n)
        Growth rates which are known. growth_data[i,0] is the density of species
        i and growth_data[i,1] is the per capita growth rate at said densities.
        Can be used to pass growth rates at to high (or low densities).
        e.g. growth_data = [[[1e20],[-dil]],
                            [[1e20],[-dil]]]
        Where dil is the dilution rate of the system.
    from_R: boolean, default False
        Set to True if function is called via R by reticulate package.
        
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
    N_t_fun: dictionary with callables
        keys are ["spec0_low", "spec1_low", "spec0_high", "spec1_high"]
        
        N_t_fun[key](t) is the fitted density of species i in the experiment 
        low respectively high. t must be a non-negative float or array
    N_t_data: array
        N_t_fun evaluated at timepoints time    
    fig: Matplotlib figure
        only returned if ``visualize`` is True
    ax: Matplotlib axes
        only returned if ``visualize`` is True
        
    Literature:
    The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
    DOI:
    """
    if from_R:
        dens = np.array([np.array(dens[0]), np.array(dens[1])])
        time = np.array(time)
        r_i = np.array(r_i)

        visualize = False # leads to fatal Error in R
    else:
        # convert input parameters
        dens = np.array(dens)
        time = np.array(time)
    input_par, log, id_exp_1 = __input_check_exp__(k , s, f0, N_star,
                    id_exp_1, log, extrapolate, dens, time, growth_data) 
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
    
    if type(N_star) == str and N_star == "spline":
        N_star = log["nanmean"](dens[...,-1], axis = (1)) # starting guess
        # compute N_star using the spline interpolations
        for i in range(2):
            N_star[i] = fsolve(f_spec, N_star[i], args = (i,))
    else:
        N_star = input_par["N_star"]
            
    # set starting guess for c
    c = np.ones((2,2))
    c[[1,0],[0,1]] = [N_star[0]/N_star[1], N_star[1]/N_star[0]] 

    # compute the ND, FD etc. parameters
    pars = {"N_star": np.array([[0,N_star[1]],[N_star[0],0]]),
            "r_i": r_i, "f": f, "c": c}
    try:
        pars = NFD_model(f, pars = pars, experimental = True)
    except InputError:
        # densities over time
        N_t_fun, N_t_data = dens_over_time(f_spec, time, dens, id_exp_1)
        fig, ax = visualize_fun(time, dens, pars, N_t_fun, 
                                    id_exp_1, log, f_spec)
        raise InputError("Fitting spline to data did not result in" 
                " reasonable results, please check input (see plots)")
        
    # densities over time
    N_t_fun, N_t_data = dens_over_time(f_spec, time, dens, id_exp_1)
    if visualize: # visualize results if necessary
        fig, ax = visualize_fun(time, dens, pars, N_t_fun, id_exp_1, log, f_spec)
        return pars, N_t_fun, N_t_data, fig, ax
    
    return pars, N_t_fun, N_t_data

def __input_check_exp__(k , s, f0, N_star, id_exp_1, log, extrapolate, 
                        dens, time, growth_data):
    # convert input parameters to right format if necessary
    
    input_par = {"k":k, "s":s}
    for key in input_par.keys():
        if type(input_par[key]) != list:
            input_par[key] = 2*[input_par[key]]
            
    input_par["extrapolate"] = extrapolate   
 
    if id_exp_1 is None: # assume first half is exp1, i.e. low starting density
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
        
    # compute data intern variance, needed for computing s
    time_diff = time[1:]-time[:-1]
    per_cap_growth = np.log(dens[...,1:]/dens[...,:-1])/time_diff
    input_par["var"] = np.nanmean(np.nanvar(per_cap_growth, axis = 1),
             axis = -1)
    
    # add predefined N_star value to growth_data if needed
    if type(N_star) == str and N_star == "average":
        N_star = log["nanmean"](dens[...,-1], axis = (1))
        if not np.all(np.isfinite(N_star)):
            raise InputError("Computed values for ``N_star`` "
                "are {} and non-finite.".format(N_star) +
                "Please pass directly or check that the last measuring time"
                "does contain at least some finite values")
    if type(N_star) != str or N_star != "spline":
        # per cap grwth is 0 at equilibrium
        growth_append = np.array([N_star, [0,0]]).T.reshape(2,2,1)
        growth_data = np.append(growth_data, growth_append, axis = -1)
    input_par["N_star"] = N_star
            
    
    # set the values for f0
    if (type(f0) == str) and f0[:5] == "perc=":
        f0 = np.nanpercentile(per_cap_growth, float(f0[5:]), axis = (1,2))
    elif type(f0) == str and f0 == "linear":
        f0 = np.nanmean(per_cap_growth[:,id_exp_1,0], axis = 1)
    if type(f0) != str and not np.all(np.isfinite(f0)):
            raise InputError("Computed values for ``f0`` "
                "are {} and non-finite.".format(f0) +
                "Please pass directly or check that the first measuring times"
                "do contain at least some finite values")
    if type(f0) != str or  f0 != "spline": # force spline to use f0
        scal = 2.0**(-np.arange(2).reshape(-1,1))
        dens_ap = scal*np.nanmin(dens, axis = (1,2))
        gr_ap = f0*np.ones((len(scal),2))
        growth_append = [[dens_ap[:,0], gr_ap[:,0]],[dens_ap[:,1], gr_ap[:,1]]]
        growth_data = np.append(growth_data, growth_append, axis = -1)        
    if np.any(np.isnan(growth_data)):
        raise InputError("Growth rates to be fitted appear to be non-finite"
            "NA are allowed in ``dens``, however not in ``time``, ``r_i``"
            ", ``f0``, ``growth_data`` or ``N_star``")        
    input_par["growth_data"] = growth_data
    return input_par, log, id_exp_1
    
    
def per_capita_growth(dens, time, input_par, log):
    """interpolate the per capita growth rate of the species
           
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
        # above maximum, use maximal entry if not extrapolate
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
    
    # use middle point for growth rate
    dens_mid = log["nanmean"]([dens[...,1:], dens[...,:-1]], axis = 0)
    
    # remove nan's
    finite = np.isfinite(dens_mid) & np.isfinite(per_cap_growth)
    dens_finite = dens_mid[finite]
    per_cap_growth = per_cap_growth[finite]
    w = np.ones(len(dens_finite)) # weighting for univariate spline
    
    # add predefined growth data
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
    if np.any(np.isnan(per_cap_growth)):
        raise InputError("Growth rates to be fitted appear to be non-finite"
            "for species {}.".format(i) +
            "NA are allowed in ``dens``, however not in ``time``, ``r_i``"
            ", ``f0`` or ``N_star``") 
    dNdt = uni_sp(log["log"](dens_finite), per_cap_growth,
                  k = input_par["k"][i], s = s, w = w)
    print(dNdt.get_coeffs(), "coeffs")
    
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
            # solve differential equation given by per-capita growth rates
            N_t_fun[key] = lambda t, i = i, N_start = N_start: odeint(
                    lambda N,t: N*f_spec(N,i), N_start,
                    np.append(time[0],t))[1:]
            N_t_data[key] = N_t_fun[key](time)[:,0]
    return N_t_fun, N_t_data
    
def visualize_fun(time, dens, pars, N_t_fun, id_exp_1, log, f_spec):
    # visualize the fitted population density and the per capita growthrate
    
    import matplotlib.pyplot as plt
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
    dens_mid = log["nanmean"]([dens[...,1:], dens[...,:-1]], axis = 0)
    for i in range(2):
        j = [1,0][i] 
        ax[1,i].plot(dens_mid[i, :, :-1], per_cap_growth[i],
          'o', color = col[i])
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
    ax[0,0].set_ylim(np.nanpercentile(dens, [0,100])*[0.8,1.2])

    return fig, ax