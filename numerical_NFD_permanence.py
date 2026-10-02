"""
@author: J.W.Spaak
Numerically compute ND and FD for a model
"""

import numpy as np
from scipy.optimize import brentq, fsolve
from scipy.integrate import solve_ivp
from warnings import warn, catch_warnings
from scipy.optimize._numdiff import approx_derivative
import matplotlib.pyplot as plt
import numbers

def NFD_model(f, n_spec = 2, args = (), monotone_f = True, pars = None,
              f_jacobian = None, xtol = 1e-5, equi_time = np.arange(1000),
              N_tol = 1e-3, f_tol = 1e-2, plot_dynamics = None, N_star = None,
              present = None, create_save_f = True, plot_results = False):
    """Compute the ND and FD for a differential equation f
    
    Compute the niche difference (ND), niche overlapp (NO), 
    fitnes difference(FD) and conversion factors (c)
    
    Parameters
    -----------
    f : callable ``f(N, *args)``
        Percapita growth rate of the species.
        1/N dN/dt = f(N)
    
    n_spec : int, optional, default = 2
        number of species in the system
    
    args : tuple, optional
        Any extra arguments to `f` and potentially to `f_jacobian`
    
    monotone_f : boolean, default = True
        Whether ``f_i(N_i,0)`` is monotonly decreasing in ``N_i``
    
    pars : dict, default {}
        A dictionary to pass arguments to help numerical solvers.
        The entries of this dictionary might be changed during the computation
    
        ``N_star`` : ndarray (shape = (n_spec, n_spec))
            N_star[i] starting guess for equilibrium density with species `i`
            absent.
        ``c`` : ndarray (shape = (n_spec, n_spec))
            Starting guess for the conversion factors from one species to the
            other. `c` is assumed to be symmetric an only the uper triangular
            values are relevant
    
    f_jacobian : callable, optional
        Jacobian of the per-capita growth function ``f``.
        Should have signature ``f_jacobian(N, *args)``.
    
    xtol : float, default 1e-5
        Precision used in solvers. This does not guarantee any final precision,
        but decreasing this value will increase precision
    
    equi_time : array-like, default = np.arange(1000)
        Time vector used for numerical integration when estimating equilibria.
        If set to a value the time will be set to no.linspace(0, equi_time, 101)
    
    N_tol : float, default 1e-3
        Minimal density below which species are assumed to be extinct in 
        community dynmics
    
    f_tol : float, default 1e-2
        Tolerance for per-capita growth rates when checking equilibrium conditions.
    
    plot_dynamics : boolean, default True
        If True, plots population dynamics during equilibrium search.
    
    N_star : ndarray, optional, lenght(n_spec)
        Initial guess for equilibrium densities of the full community
        If None, equilibrium densities are estimated numerically.
    
    present : boolean array, optional
        Boolean array (length n_spec) indicating which species are initially present.
        If None, all species are assumed present for the computation of niche
        and fitness differences
    
    create_save_f : boolean, default True
        Creates a wrapper around the growth rate that checks for non-finite values.
        Setting this to False makes the computation faster,
        but makes debugging harder.
    
    plot_results : boolean, default False
        Set to True if you want a plot of the final niche and fitness differences.
        
    Returns
    -------
    pars : dict
        A dictionary with the following keys: 
            
    ``ND`` : ndarray (shape = n_spec)
        Niche difference of the species to the other species
        ND = (r_i - eta)/(\mu -eta)
    ``FD`` : ndarray (shape = n_spec)
        Fitness difference according to Spaak and De Laender 2020
        FD = fc/f0
    ``c`` : ndarray (shape = (n_spec, n_spec))
        The conversion factors from one species to the other.    
    ``N_star`` : ndarray (shape = (n_spec, n_spec))
        N_star[i] equilibrium density with species `i` absent. N_star[i,i] is 0
    ``r_i`` : ndarray (shape = n_spec)
        invasion growth rates of the species
    ``eta``: ndarray (shape = n_spec)
        no-niche growth rate f(\sum c_j^i N_j^(-i),0)
    ``mu_i``: ndarray (shape = n_spec)
        intrinsic growth rate f(0,0)
    ``surv``: ndarray (shape = n_spec, dtype = bool)
        Array indicating which species survives
    ``present``: ndarray (shape = n_spec, dtype = bool)
        Which species are assumed to be present in the community for computation
        Equal to input parameter ``present``
    ``f``: callable
        per-capita growth rate, passed for convenience
    ``f_jac``: callable
        jacobian of per-capita growth rate, passed for convenience
        ``f_jac`` is equal to ``f_jacobian``
    
    Raises:
        InputError:
            Is raised if system cannot automatically solve equations.
            Starting estimates for N_star and c should be passed.
    
    Examples:
        See "Example,compute NFD.py" and "Complicated examples for NFD.py"
        for applications for models
        See "Exp_plots.py" for application to experimental data
    
    Debugging:
        If InputError is raised the problem causing information is saved in
        pars.
        To access it rerun the code in the following way (or similar)
            
        pars = {}
        pars = NFD_model(f, pars = pars)
        print(pars)
        
        pars will then contain additional information  
        
    Literature:
    "Intuitive and broadly applicable definitions of 
    niche and fitness differences", J.W.Spaak, F. deLaender
    DOI: https://doi.org/10.1101/482703 
    """

    # check input on correctness
    (f, n_spec, args, monotone_f, pars, f_jacobian, xtol,
     equi_time, N_tol, f_tol, plot_dynamics, N_star,
     present) = __input_check__(f, n_spec, args, monotone_f, pars, f_jacobian, xtol,
      equi_time, N_tol, f_tol, plot_dynamics, N_star, present)
    
    equi_args = {"time": equi_time,
                 "f_tol": f_tol,
                 "N_tol": N_tol,
                 "xtol": xtol,
                 "plot_dynamics": plot_dynamics,
                 "N_star": N_star}
    
    # obtain equilibria densities and invasion growth rates    
    pars = preconditioner(f, f_jacobian, args, n_spec, pars, xtol, monotone_f,
                          create_save_f, present)    
             
    # find equilibrium community
    equi = find_equi(equi_args["N_star"][pars["present"]], pars["present"],
                     pars, equi_args)
    
    # identify which species survived
    pars["surv"] = equi["equi"] != 0
    pars["equi"] = equi["equi"]
    
    if not any(pars["surv"]):
        raise InputError("Apparently no species survives in the full community."
              " Please check your model. Potentially rerun the `NFD_model`"
              " with `plot_dynamics = True` for more info.")
        
    # add resident community for all non-surviving species
    pars["N_star"][~pars["surv"]] = np.tile(pars["equi"], (sum(~pars["surv"]),1))
    
    # indices of all surviving species
    ind_surv = np.where(pars["surv"])[0]
    # find resident communities for surviving species
    for i in ind_surv:
        pres_temp = pars["surv"].copy()
        pres_temp[i] = False
        if any(pres_temp): # more than one species survive
            equi_temp = find_equi(pars["N_star"][i, pres_temp],
                                  pres_temp, pars, equi_args)
            pars["N_star"][i] = equi_temp["equi"]
        else:
            pars["N_star"][i] = 0
    
    # check whether all densities are above threshhold
    if np.any((pars["N_star"]>0) & (pars["N_star"]<N_tol)):
           N_tol_suggest = np.amin(pars["N_star"][pars["N_star"]>0])*0.9
           warn(f"Some species have equilibrium densities below proposed "
                f"minimal density `N_tol` = {N_tol}. This may be fine, but "
                f"please double check. Decreaseing `N_tol` to {N_tol_suggest}"
                f" will likely get rid of this warning.")
    
    # compute the invasion growth rates for each species
    pars["r_i"] = np.array([pars["f"](pars["N_star"][i])[i] for i in range(n_spec)])
    # check that invasion matches coexistence requirements
    if any((pars["r_i"][present]>0) != pars["surv"][present]):
        raise InputError("Automatically computed equilibria lead to wrong coexistence. "
                         "To solve this problem please provide equilibrium estimates. "
                         "via `pars['N_star'] argument. Alternatively, please"
                         "increase `equi_time`, decrease `N_tol` and/or "
                         "decrease `f_tol`.")
    
    # list of all species
    l_spec = list(range(n_spec))
    # compute conversion factors
    for i in l_spec:
        for j in l_spec:
            if i>=j: # c is assumed to be symmetric, c[i,i] = 1
                continue
            if pars["N_star"][i,j] == pars["N_star"][j,i] == 0:
                pars["c"][[i,j],[j,i]] = np.nan
            else:
                pars["c"][[i,j],[j,i]] = solve_c(pars,[i,j],monotone_f,xtol=xtol)



    # compute no-niche growth rates
    pars["eta_i"] = np.empty(n_spec)
    # compute intrinsic growth rates
    pars["mu_i"] = pars["f"](np.zeros(n_spec))
    
    for i in range(n_spec):
        # creat a list with i at the beginning [i,0,1,...,i-1,i+1,...,n_spec-1]
        sp = np.array([i]+l_spec[:i]+l_spec[i+1:])
        
        pars["eta_i"][i] = pars["f"](switch_niche(pars["N_star"][i],sp,pars["c"][i, sp[1:]]))[i]

    with catch_warnings(record = True):
        pars["ND"] = (pars["r_i"] - pars["eta_i"])/(pars["mu_i"] - pars["eta_i"])
        pars["FD"] = (0 - pars["eta_i"])/(pars["mu_i"] - pars["eta_i"])
        
        # for species with mu = eta = r_i
        pars["FD"][np.isnan(pars["ND"])] = 1
        pars["ND"][np.isnan(pars["ND"])] = 1
    
    if plot_results:
        plot_NFD(pars)

    return pars  

def find_equi(N_start, present, pars, equi_args):
    
    # define which community we're solving
    com = ("full-community" if all(present == pars["present"])
           else "-{} community".format(np.where(present != pars["surv"])[0][0]))
    
    # make sure every species that is present has positive starting density
    if np.any(N_start<=0):
        N_start_replace = 1 if np.all(N_start<=0) else np.mean(N_start[N_start>0])
        warn(f"In {com} species {np.arange(len(present))[present][N_start<=0]} had non-positive"
             "starting densities. These have been overwritten automatically to "
             f"{N_start_replace}")
        N_start[N_start<=0] = N_start_replace
        
        
    # vector containing all species
    N_all = np.zeros(pars["n_spec"])
    
    hill_exp = 4
    
    # convert per-capita growth rate to actual growth rate
    def ode_fun(t, N):
        # add densities of absent species
        N_all[present] = N
        return pars["f"](N_all)[present] * N      # only present species
    
    # event function for solve_ivp (stop when stop_ode == 0)
    def event(t, N):
        # add densities of absent species
        N_all[present] = N
        
        # do not focus on extinct species
        weight = N_all**hill_exp/(N_all**hill_exp + (equi_args["N_tol"]/2)**hill_exp)
        return np.max(np.abs(weight*pars["f"](N_all))) - equi_args["f_tol"]  
    
    event.terminal = True
    event.direction = -1
    
    
    if event(0, N_start) <0:
        # initial estimate already close to equilibrium
        N_all[present] = N_start
        sol = None # define sol for consistency of later code
    else:
        # solve differential equation to equilibrium
        sol = solve_ivp(
            ode_fun,
            t_span=(equi_args["time"][0], equi_args["time"][-1]),
            y0=N_start, dense_output=True,
            events=event, atol = equi_args["xtol"], method = "LSODA"
        )
        
        
        
        if len(sol.t_events[0]) > 0:
            # add the potential stopping event to the densities over time
            sol.t = np.append(sol.t, sol.t_events[0])
            sol.y = np.append(sol.y, sol.y_events[0].T, axis = 1)
        
        # add the densities at the equi_args["time"]
        sol.t = np.append(sol.t, equi_args["time"][equi_args["time"] < max(sol.t)])
        ind = np.argsort(sol.t)
        sol.y = np.append(sol.y, sol.sol(equi_args["time"]), axis = 1)
        sol.y = sol.y[:,ind]
        sol.t = sol.t[ind]

        N_all[present] = sol.y[:, -1]
        

        if not sol.success and (not equi_args["plot_dynamics"] == False):
            plot_equi(sol, N_all[present], present[present], equi_args, com)
            raise RuntimeError(f"Numerical integration failed for community {com}")
  
    # remove all absent species from equilibrium
    present_new = N_all>equi_args["N_tol"]
    N_all[~present_new] = 0 # remove species that went extinct
    # did at least one species survive?
    if not np.any(present_new):
    
        warn(f"All species went extinct for {com}.")
    
        if ((equi_args["plot_dynamics"] or # plot forced
            np.any(pars["f"](N_all)[present] > equi_args["xtol"])) # error present
            and (not equi_args["plot_dynamics"] == False)): # plotting not turned off
            plot_equi(sol, N_all[present], present_new[present], equi_args, com)              
    
        if np.any(pars["f"](N_all)[present] > equi_args["xtol"]):
             raise RuntimeError(
                f"Could not find automatical equilibrium for {com}. "
                "Please check plot. To solve this problem, please provide "
                "equilibrium estimates via `pars['N_star']` or `N_star` argument. "
                "Alternatively, please increase equi_time, decrease N_min, "
                "and/or decrease dN_dt."
            )
    
        return {"equi": N_all, "equi_found": True}
    
    
    # --- root finder refinement ---
    def fsolve_fun(N):
        N_all[present_new] = N
        # use actual growth rates, not per-capita
        return(N*pars["f"](N_all)[present_new]) # only return present species
    
    try:
        x, info, success, msg = fsolve(fsolve_fun, N_all[present_new],
                                       xtol=equi_args["xtol"],
                                       full_output = True)
        if success == 1:
            # update equilibrium density
            N_all[present_new] = x
            
            # Check stability of equilibrium
            # Jacobian of system at equilibrium
            r = np.zeros((sum(present_new), sum(present_new)))
            r[np.triu_indices(sum(present_new))] = info["r"].copy()
            jac = info["fjac"].T.dot(r)
        else:
            jac = approx_derivative(fsolve_fun, N_all[present_new])
            
    except:
        # can't update equilibrium
        # compute jacobian
        jac = approx_derivative(fsolve_fun, N_all[present_new])  
    
    # species that were present but went extinct
    extinct = ~present_new.copy()
    extinct[~present] = False
    
    # check whether equilibrium has been reached
    f_vals = pars["f"](N_all)
    
    equi_found = (
        np.max(np.abs(f_vals[present_new]) < equi_args["xtol"])  # equilibrium
        and np.all(N_all[present_new] > 0)                       # feasible
        and np.all(np.isfinite(N_all))                           # valid numbers
        and np.max(np.real(np.linalg.eigvals(jac))) < 0          # stable
        and np.all(f_vals[extinct] < 0)                          # no invasion
    )
    
    
    if (equi_args["plot_dynamics"] or (not equi_found)) and (not equi_args["plot_dynamics"] == False):
        plot_equi(sol, N_all[present], present_new[present], equi_args, com)
    
    
    if not equi_found:
        pars["sol"] = sol
        pars["com"] = com
        pars["temp_equi"] = N_all
        pars["present"] = present
        raise InputError(
            f"Could not find automatical equilibrium for {com}. "
            "Please check plot. To solve this problem, please provide "
            "equilibrium estimates via `pars['N_star']` or `N_star` argument. "
            "Alternatively, please increase equi_time, decrease N_tol, "
            "and/or decrease f_tol."
        )
    
    
    return {
        "equi": N_all,
        "equi_found": equi_found,
        "jac": jac,
        "no_invasion": f_vals,
    }

def plot_equi(sol, equi, loc_surv, equi_args, com):
    if sol is None:
        warn(
            f"Can't plot dynamics because species already at equilibrium for {com}."
        ) # sol was never computed, no plots possible
        return
    # layout
    plt.figure()
    plt.xlim(sol.t[[0,-1]])
    plt.ylim(equi_args["N_tol"] / 100, np.max(sol.y) * 1.5)
    plt.yscale("log")

    plt.title(f"{com}")
    plt.xlabel("Time")
    plt.ylabel("Density")

    # plot densities over time
    for i in range(sol.y.shape[0]):
        color = "black" if loc_surv[i] else "red"
        plt.plot(sol.t, sol.y[i], color=color)
        plt.scatter(sol.t[-1], equi[i], color=color)

    plt.axhline(equi_args["N_tol"], color="red")
        
class InputError(Exception):
    pass
        
def preconditioner(f, f_jacobian, args, n_spec, pars, xtol, monotone_f,
                   create_save_f, present):
    """Returns equilibria densities and invasion growth rates for system `f`
    
    Parameters
    -----------
    same as `find_NFD`
            
    Returns
    -------
    pars : dict
        A dictionary with the keys:
        
        ``N_star`` : ndarray (shape = (n_spec, n_spec))
            N_star[i] is the equilibrium density of the system with species 
            i absent. The density of species i is set to 0.
        ``r_i`` : ndarray (shape = n_spec)
            invsaion growth rates of the species
    """        
    # expected shapes of pars
    pars_def = {"c": np.ones((n_spec,n_spec)),
                "r_i": np.zeros(n_spec),
                "N_star": np.ones((n_spec, n_spec))}
    
    warn_string = "pars[{}] must be array with shape {}."\
                +" The values will be computed automatically"
    # check given keys of pars for correctness
    for key in pars_def.keys():
        try:
            if pars[key].shape == pars_def[key].shape:
                pass
            else: # `pars` doesn't have expected shape
                pars[key] = pars_def[key]
                warn(warn_string.format(key,pars_def[key].shape))
        except KeyError: # key not present in `pars`
            pars[key] = pars_def[key]
        except AttributeError: #`pars` isn't an array
            pars[key] = pars_def[key]
            warn(warn_string.format(key,pars_def[key].shape))         
    
    # create a save function f
    if create_save_f:
        def save_f(N):
            N_save = N.copy()
            N_save[N_save<0] = 0
            
            # negative infinite growth assumed at infinite densities
            if np.any(~np.isfinite(N)):
                return np.full(n_spec, -np.inf)
            
            try:
                growth = f(N_save, *args)
                if np.any(~np.isfinite(growth)):
                    raise InputError(f"Function call resulted in non-finite growth at densities {N}")
            except:
                warn(f"Function call of f did not work properly at densities {N}")
                raise
            return growth
        pars["f"] = save_f
        
        if f_jacobian is None:
            pars["f_jac"] = None
        else:
            def save_f_jacobian(N):
                try:
                    jacobian = f(N, *args)
                    if np.any(~np.isfinite(jacobian)):
                        raise InputError(f"Function call of jacobian resulted in non-finite values at densities {N}")
                except:
                    warn(f"Function call of jacobian did not work properly at densities {N}")
                    raise
                return jacobian
            pars["f_jac"] = None
        
    else:
        pars["f"] = lambda N: f(N, *args)
        # test function once
        try:
            N = np.random.uniform(0,1, n_spec)
            pars["f"](N)
        except:
            warn(f"Function call of f did not work properly at densities {N}")
            raise
        
        if f_jacobian is None:
            pars["f_jac"] = None
        else:
            # test jacobian once
            pars["f_jac"] = lambda N: f_jacobian(N, *args)
            try:
                pars["f_jac"](N)
            except:
                warn(f"Function call of jacobian did not work properly at densities {N}")
                raise
        
            
    # c must be a positive real number
    if (np.any(~np.isfinite(pars["c"])) or np.any(pars["c"]<=0) 
            or pars["c"].dtype != float):
        warn("Some entries in pars['c'] were not positive real numbers."
             "These are replaced with 1")
        pars["c"] = np.real(pars["c"])
        pars["c"][pars["c"] <= 0] = 1
        pars["c"][~np.isfinite(pars["c"])] = 1
    
    pars["present"] = present
    pars["n_spec"] = n_spec

    return pars
    
def solve_c(pars, sp = [0,1], monotone_f = True, xtol = 1e-10):
    """find the conversion factor c for species sp
    
    Parameters
    ----------
    pars : dict
        Containing the N_star and r_i values, see `preconditioner`
    sp: array-like
        The two species to convert into each other
        
    Returns
    -------
    c : float, the conversion factor c_sp[0]^sp[1]
    """
    NO_funs = [NO_fun(pars, sp),
               NO_fun(pars, sp[::-1])]
    NO_values = [NO_funs[0](pars["c"][sp[0], sp[1]]),
                 NO_funs[1](pars["c"][sp[1], sp[0]])]
    # do species interact?
    if np.isclose(NO_values, [0,0]).any():
        return special_case(np.isclose(NO_values, [0,0]), sp)
    
    # has one species reached minimal growth rate?
    if np.isinf(NO_values).any():
        return special_case_mort(np.isinf(NO_values), sp)
    
    sp = np.asarray(sp)
    
    def inter_fun(c):
        # equation to be solved
        NO_ij = np.abs(NO_funs[0](c))
        NO_ji = np.abs(NO_funs[1](1/c))
        return NO_ij-NO_ji
    
    # use a generic numerical solver when `f` is not montone
    # potentially there are multiple solutions
    if not monotone_f:
        c = fsolve(inter_fun,pars["c"][sp[0],sp[1]],xtol = xtol)[0]
        if np.abs(inter_fun(c))>xtol:
            pars["c found by fsolve"] = c
            raise InputError("Not able to find c_{}^{}.".format(*sp) +
                "Please pass a better guess for c_i^j via the `pars` argument")
        return c, 1/c
        
    # if `f` is monotone then the solution is unique, find it with a more
    # robust method
        
    # find interval for brentq method
    a = pars["c"][sp[0],sp[1]]
    # find which species has higher NO for c0
    direction = np.sign(inter_fun(a))
    
    if direction == 0: # starting guess for c is correct
        return a, 1/a
    fac = 2**direction
    if not np.isfinite(direction):
        pars["function inputs"] = [switch_niche(pars["N_star"][es[0]],es,c)
                for c in [0,a, 1/a] for es in [sp, sp[::-1]]]
        pars["function outputs"] = [pars["f"](inp) for 
             inp in pars["function inputs"]]
        raise InputError("function `f` seems to be returning nonfinite values")
    b = float(a*fac)
    # change searching range to find c with changed size of NO
    while np.sign(inter_fun(b)) == direction:
        a = b
        b *= fac
        # test whether a and be behave as they should (e.g. nonfinite)
        if not((2*a == b) or (2*b == a)) or np.sign(b-a) != direction:
            raise InputError("Not able to find c_{}^{}.".format(*sp) +
                "Please pass a better guess for c_i^j via the `pars` argument"+
                ". Please also check for non-positive entries in pars[``c``]")
    # solve equation
    try:
        c = brentq(inter_fun,a,b)
    except ValueError:
        raise ValueError("f does not seem to be monotone. Please run with"
                         +"`monotone_f = False`")
    # test whether c actually is correct
    # c = 0 implies issue with brentq
    if (c==0) or inter_fun(c)>xtol:
        pars["c"][sp[0],sp[1]] = c
        raise InputError("Not able to find c_{}^{}.".format(*sp) +
                "Please pass a better guess for c_i^j via the `pars` argument"+
                ". Please also check for non-positive entries in pars[``c``]")
    return c, 1/c # return c_i and c_j = 1/c_i

def NO_fun(pars, sp):
    # return a function that computes the niche overlap as a function of c
    
    if pars["N_star"][sp[0], sp[1]] > 0: # competitor present
        f0 = pars["f"](switch_niche(pars["N_star"][sp[0]],sp))[sp[0]]
        def NO_pre(c_ij):
            fc = pars["f"](switch_niche(pars["N_star"][sp[0]],sp,c_ij))[sp[0]]
            return np.abs( (f0 - pars["r_i"][sp[0]])/(f0-fc))
    else: # competitor absent, use jacobian
        if pars["f_jac"] is None: # compute jacobian numerically
            jac = approx_derivative(pars["f"], pars["N_star"][sp[0]])
        else:
            jac = pars["f_jac"](pars["N_star"][sp[0]])
        def NO_pre(c_ij):
            return 1/c_ij*jac[sp[0], sp[1]]/jac[sp[0], sp[0]]
    
    return lambda c_ij: NO_pre(c_ij) if np.isfinite(NO_pre(c_ij)) else 1

def special_case(no_comp, sp):
    # Return c for special case where one spec is not affected by competition
    
    warn("Species {} and {} do not seem to interact.".format(sp[0], sp[1]) +
      " This may result in nonfinite c, ND and FD values.")
    
    if no_comp.all():
        return 0, 0 # species do not interact at all, c set to zero
    elif (no_comp == [True, False]).all():
        return 0, np.inf # only first species affected
    elif (no_comp == [False, True]).all():
        return np.inf, 0
    
def special_case_mort(mort, sp):
    # Return c for special case where one spec is not affected by itself
    
    warn("Species {} or {} reached mortality rate.".format(sp[0], sp[1]) +
      " This may result in nonfinite c, ND and FD values.")
    
    if mort.all():
        return 0, 0 # both species have reached mortality rate
    elif (mort == [True, False]).all():
        return np.inf, 0 # only first species affected
    elif (mort == [False, True]).all():
        return 0, np.inf
    
def switch_niche(N,sp,c=0):
    # switch the niche of sp[1:] into niche of sp[0]
    N = N.copy()
    N[sp[0]] += np.nansum(c*N[sp[1:]])
    N[sp[1:]] = 0
    return N

def __input_check__(
    f,
    n_spec=2,
    args=(),
    monotone_f=True,
    pars=None,
    f_jacobian=None,
    xtol=1e-5,
    equi_time=None,
    N_tol=1e-3,
    f_tol=1e-2,
    plot_dynamics=True,
    N_star=None,
    present=None, create_save_f = True
    ):


    # --- f ---
    if not callable(f):
        raise TypeError("f must be callable")
    
    # --- n_spec ---
    if not (isinstance(n_spec, int) and n_spec > 0):
        raise ValueError("n_spec must be a positive integer")

    # --- args ---
    if not isinstance(args, tuple):
        raise TypeError("args must be a tuple")
      
    # try whether function call actually works
    try:
        f(np.ones(n_spec), *args)
    except:
        warn("Function call did not work properly")
        raise

    # --- monotone_f ---
    if not isinstance(monotone_f, bool):
        raise TypeError("monotone_f must be boolean")

    # --- pars ---
    if pars is None:
        pars = {}
    if not isinstance(pars, dict):
        raise TypeError("pars must be a dictionary")

    # --- f_jacobian ---
    if f_jacobian is None:
        pass
    elif not callable(f_jacobian):
        raise TypeError("f_jacobian must be callable or None")
    else:
        # try function call of jacobian
        try:
            f_jacobian(np.ones(n_spec), *args)
        except:
            warn("Function call did not work properly")
            raise

    # --- xtol ---
    if not (isinstance(xtol, numbers.Real) and xtol > 0):
        raise ValueError("xtol must be positive")

    # --- equi_time ---
    if equi_time is None:
        equi_time = np.linspace(0, 1000, 101)

    if isinstance(equi_time, numbers.Real):
        if equi_time <= 0:
            raise ValueError("equi_time must be positive")
        equi_time = np.linspace(0, equi_time, 101)
    else:
        equi_time = np.asarray(equi_time)
        if (
            equi_time.ndim != 1
            or not np.all(np.diff(equi_time) > 0)
        ):
            raise ValueError(
                "equi_time must be increasing positive numbers"
            )

    # --- N_tol ---
    if not (isinstance(N_tol, numbers.Real) and N_tol > 0):
        raise ValueError("N_tol must be positive")

    # --- f_tol ---
    if not (isinstance(f_tol, numbers.Real) and f_tol > 0):
        raise ValueError("f_tol must be positive")

    # --- plot_dynamics ---
    if not (isinstance(plot_dynamics, bool) or (plot_dynamics is None)):
        raise TypeError("plot_dynamics must be boolean or None")

    # --- N_star ---
    if N_star is None:
        N_star = np.ones(n_spec)
    else:
        N_star = np.asarray(N_star)
        if N_star.shape != (n_spec,):
            raise ValueError("N_star must be length n_spec")

    # --- present ---
    if present is None:
        present = np.ones(n_spec, dtype=bool)
    else:
        present = np.asarray(present, dtype=bool)
        if present.shape != (n_spec,):
            raise ValueError("present must be boolean array of length n_spec")

    return (
        f,
        n_spec,
        args,
        monotone_f,
        pars,
        f_jacobian,
        xtol,
        equi_time,
        N_tol,
        f_tol,
        plot_dynamics,
        N_star,
        present,
    )

def plot_NFD(pars, ax = None):
    if ax is None:
        fig = plt.figure()
        ax = plt.gca()
    

    try:
        names = pars["names"]
    except KeyError:
        names = np.arange(pars["n_spec"])
        
    colors = np.where(pars["surv"], "g", "r")
    colors[~pars["present"]] = "b"
    for i in range(pars["n_spec"]):
        ax.text(pars["ND"][i], pars["FD"][i], names[i],
                color = colors[i], ha = "center", va = "center")
    
    # add layout
    ax.set_xlabel("Niche differences")
    ax.set_ylabel("Fitness differences")
    x_range = [min(pars["ND"][np.isfinite(pars["ND"])]),
               max(pars["ND"][np.isfinite(pars["ND"])])]
    ax.set_xlim(np.mean(x_range) + np.array([-0.6, 0.6])*(x_range[1]-x_range[0]))
    y_range = [min(pars["FD"][np.isfinite(pars["FD"])]),
               max(pars["FD"][np.isfinite(pars["FD"])])]
    ax.set_ylim(np.mean(y_range) + np.array([-0.6, 0.6])*(y_range[1]-y_range[0]))
    
    labels = {"r": "extinct", "g": "Persists", "b": "Not in community"}
    # add legend
    for col in np.unique(colors):
        ax.plot(np.nan, np.nan, 'o', color = col,
                label = labels[col])
    ax.legend()
    
    # add boundaries
    ax.axvline(1, color = "grey")
    ax.axvline(0, color = "grey")
    ax.axhline(0, color = "grey")
    ax.axhline(1, color = "grey")
    ax.plot(ax.get_xlim(), ax.get_xlim(), 'k')
    
    return fig, ax


if False and __name__ == "__main__":
    n_spec = 20
    # create random interaction matrix
    A = np.random.uniform(-0.2,0.6, (n_spec, n_spec))
    np.fill_diagonal(A, np.random.uniform(0.9, 1.2, n_spec))

    mu = np.random.uniform(0.9,1.1, n_spec)
    f = lambda N: 2- A.dot(N)
    # define all variables
    n_spec = len(A) 
    args = ()
    monotone_f = True
    pars = None
    f_jacobian = None
    xtol = 1e-5
    equi_time = np.arange(1000)
    N_tol = 1e-3
    f_tol = 1e-2
    plot_dynamics = True
    N_star = None
    present = None
    create_save_f = True
    
    def LV_model(N, mu, A):
        return (mu - A.dot(N))
    pars = NFD_model(LV_model, n_spec = len(A), args = (mu, A),
                     f_tol = 1e-4, plot_results=True)
    
    ###########################################################################
    # check results
    r_i_all = np.empty((n_spec, n_spec))
    for i in range(n_spec):
        r_i_all[i] = LV_model(pars["N_star"][i], mu, A)
    
    # growth rate of present species must be small
    assert(np.amax(np.abs(r_i_all[pars["N_star"]>0]))<1e-5)
    
    # extinct species must have negative growth rates
    extinct = pars["N_star"] == 0 # extinct species
    np.fill_diagonal(extinct, False) # species -i communities are not considered extinct
    extinct[:, ~pars["surv"]] = False # died out species to not count
    
    assert(np.all(r_i_all[extinct] < 0))

    
    # invasion growth rates must match the survival pattern
    assert(np.all(np.diag(r_i_all>0) == pars["surv"]))
    
    # compute conversion factors
    c_ij = np.sqrt(np.abs(A/A.T/np.diag(A)[:, np.newaxis]*np.diag(A)))
    loc = np.where(np.abs(c_ij - pars["c"])>1e-2)

    
    assert(np.allclose(c_ij[np.isfinite(pars["c"])],
                       pars["c"][np.isfinite(pars["c"])]))
    
    # compute eta manually
    eta_i = mu - np.diag(A)*np.einsum("ij, ij->i", pars["N_star"], c_ij)
    assert(np.allclose(eta_i, pars["eta_i"]))
    # check that invasion growth rate is computed correctly
    assert(np.allclose(np.diag(r_i_all), pars["r_i"]))
    # check intrinsic growth rate
    assert(np.allclose(mu, pars["mu_i"]))
    
    # compute niche and fitness differences
    ND = (np.diag(r_i_all) - eta_i)/(mu - eta_i)
    FD = (0 - eta_i)/(mu - eta_i)
    
    assert(np.allclose(ND, pars["ND"]))
    assert(np.allclose(FD, pars["FD"]))
    
    print("Niche and fitness differences computed correctly")
    