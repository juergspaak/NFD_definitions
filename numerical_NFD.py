"""
@author: J.W.Spaak
Numerically compute ND and FD for a model
"""

import numpy as np
from scipy.optimize import brentq, fsolve
from warnings import warn

def NFD_model(f, n_spec = 2, args = (), monotone_f = True, pars = None,
             experimental = False, from_R = False, xtol = 1e-5,
             estimate_N_star_mono = False):
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
        Any extra arguments to `f`
    monotone_f : boolean or array of booleans (lenght: n_spec), default = True
        Whether ``f_i(N_i,0)`` is monotonly decreasing in ``N_i``
        Can be specified for each function separatly by passing an array.
    pars : dict, default {}
        A dictionary to pass arguments to help numerical solvers.
        The entries of this dictionary might be changed during the computation
        
        ``N_star`` : ndarray (shape = (n_spec, n_spec))
            N_star[i] starting guess for equilibrium density with species `i`
            absent. N_star[i,i] is set to 0 
        ``r_i`` : ndarray (shape = n_spec)
            invsaion growth rates of the species
        ``c`` : ndarray (shape = (n_spec, n_spec))
            Starting guess for the conversion factors from one species to the
            other. `c` is assumed to be symmetric an only the uper triangular
            values are relevant
    experimental: boolean, default False
        Automatically set to True when used in combination with data of
        experiments. Do not set this to True manually!        
    from_R: boolean, default False
        Set to True if function is called via R by reticulate package.
        Converts types of f and equilibria.
    xtol: float, default 1e-10
        Precision requirement of solving
    estimate_N_star_mono: boolean, default False
        If True, then N_star[i,j] will be estimated with monoculture 
        equilibrium density of species j.
        Setting to True will potentially reduce speed, but result in more
        robust behaviour.
        Can only be used if ``f`` is monotone, i.e. monotone_f == True
        
    Returns
    -------
    pars : dict
        A dictionary with the following keys: 
            
    ``N_star`` : ndarray (shape = (n_spec, n_spec))
        N_star[i] equilibrium density with species `i`
        absent. N_star[i,i] is 0
    ``r_i`` : ndarray (shape = n_spec)
        invasion growth rates of the species
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
    ``fc``: ndarray (shape = n_spec)
        no-niche growth rate f(\sum c_j^i N_j^(-i),0)
    
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
    if n_spec == 1:
        # species case, single species are assumed to have ND = 1
        raise InputError("ND and FD are not (properly) defined for a single"
            "species community."
            "If needed assign manualy ND = 1 and FD = 0 for this case")
    
    if from_R:
        if n_spec-int(n_spec) == 0:
            n_spec = int(n_spec)
        else:
            raise InputError("Number of species (`n_spec`) must be an integer")
        fold = f
        #f(0)
        def f(N, *args):
            # translate dataframes, matrices etc to np.array
            return np.array(fold(N, *args)).reshape(-1)
        
        if not(pars is None):
            try:
                for key in pars.keys(): # convert to np array and make writable
                    pars[key] = np.array(pars[key])
            except AttributeError:
                raise InputError("Argument ``pars`` must be a dictionary or a"
                    "labeled list. e.g. ``pars = list(N_star = N_star)")
    # check input on correctness
    monotone_f = __input_check__(n_spec, f, args, monotone_f, pars)
    
    if experimental:
        if not ("c" in pars.keys()):
            pars["c"] = np.ones((n_spec, n_spec))
        if not ("r_i" in pars.keys()):
            pars["r_i"] = np.array([f(pars["N_star"][i], *args)[i] 
                        for i in range(n_spec)])
    if not experimental:
        # obtain equilibria densities and invasion growth rates    
        pars = preconditioner(f, args,n_spec, pars, xtol, monotone_f,
                              estimate_N_star_mono)                  
    # list of all species
    l_spec = list(range(n_spec))
    # compute conversion factors
    c = np.ones((n_spec,n_spec))
    for i in l_spec:
        for j in l_spec:
            if i>=j: # c is assumed to be symmetric, c[i,i] = 1
                continue
            c[[i,j],[j,i]] = solve_c(pars,[i,j],
                         monotone_f[i] and monotone_f[j],xtol=xtol)

    # compute NO and FD
    NO = np.empty(n_spec)
    FD = np.empty(n_spec)
    
    for i in l_spec:
        # creat a list with i at the beginning [i,0,1,...,i-1,i+1,...,n_spec-1]
        sp = np.array([i]+l_spec[:i]+l_spec[i+1:])
        # compute NO and FD
        if (c[i, sp[1:]] == 0).all():
            NO[i] = 0 # species does not interact with each other species
        else:
            NO[i] = NO_fun(pars, c[i, sp[1:]], sp)
        FD[i] = FD_fun(pars, c[i, sp[1:]], sp)
    
    # prepare returning values
    pars["NO"] = NO
    pars["ND"] = 1-NO
    pars["FD"] = FD
    pars["c"] = c
    pars["f0"] = pars["f"](np.zeros(n_spec)) # monoculture growth rate
    pars["fc"] = FD*pars["f0"] # no niche growth rate
    return pars
  
def __input_check__(n_spec, f, args, monotone_f, pars):
    # check input on (semantical) correctness
    if not isinstance(n_spec, int):
        raise InputError("Number of species (`n_spec`) must be an integer")
    
    # check whether `f` is a function and all species survive in monoculutre
    try:
        f0 = f(np.zeros(n_spec), *args)
        if f0.shape != (n_spec,):
            if not (pars is None):
                pars["function_call"] = "f(0)"
                pars["return_value"] = f0
            raise InputError("`f` must return an array of length `n_spec`")   
    except TypeError:
        print("function call of `f` did not work properly")
        raise
    except AttributeError:
        fold = f
        f = lambda N, *args: np.array(fold(N, *args))
        f0 = f(np.zeros(n_spec), *args)
        warn("`f` does not return a proper `np.ndarray`")

    # broadcast monotone_f if necessary
    return np.logical_and(monotone_f, np.full(n_spec, True, bool))
        
class InputError(Exception):
    pass
        
def preconditioner(f, args, n_spec, pars, xtol, monotone_f,
                   estimate_N_star_mono):
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
    if pars is None:
        pars = {}
        
    # expected shapes of pars
    pars_def = {"c": np.ones((n_spec,n_spec)),
                "r_i": np.zeros(n_spec)}
    
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
    
    # estimate N_star as monoculture equilibrium densities
    try:
        if pars["N_star"].shape == (n_spec, n_spec):
            pass # correct shape
        elif pars["N_star"].shape == (n_spec):
            # assume pars["N_star"][i] is monoculture of species i
            pars["N_star"] = pars["N_star"]*np.ones((n_spec, n_spec))
        else: # `pars` doesn't have expected shape
            pars["N_star"] = np.ones(n_spec)
            estimate_N_star_mono = True
            warn(warn_string.format(key,pars_def[key].shape))
    except KeyError: # key not present in `pars`
        pars["N_star"] = np.ones(n_spec)
        estimate_N_star_mono = True
    except AttributeError: #`pars` isn't an array
        pars[key] = pars_def[key]
        warn(warn_string.format(key,pars_def[key].shape))            
            
    def save_f(N):
            # allow passing infinite species densities to per capita growthrate
            if np.isinf(N).any():
                return np.full(N.shape, -np.inf)
            else:
                N = N.copy()
                N[N<0] = 0 # function might be undefined for negative densities
                return f(N, *args)
    pars["f"] = save_f
    # monoculture growth rate
    pars["f0"] = pars["f"](np.zeros(n_spec))
    
    if estimate_N_star_mono:
        if np.ndim(pars["N_star"])==2:
            N_star_mono = np.mean(pars["N_star"], axis = 0)
        else:
            N_star_mono = pars["N_star"]
        
        # starting estimates for N_star must be positive real numbers
        N_star_mono[N_star_mono<0] = 1
        N_star_mono[~np.isfinite(N_star_mono)] = 1
        
        # estimate N_star via brentq algorithm
        for i in range(n_spec):
            if pars["f0"][i]<0:
                continue # species can't survive in monoculture
                
            if monotone_f[i]:
                
                counter = 0
                # to avoid overflow, how often can we double?
                max_counter = (np.log(np.finfo(float).max) - 
                               np.log(N_star_mono[i]))/np.log(2)
                growth = pars["f"](np.insert(np.zeros(n_spec-1), i,
                                 N_star_mono[i]))[i]
                while (growth>0 and counter < max_counter-2):
                    N_star_mono[i] *= 2
                    growth = pars["f"](np.insert(np.zeros(n_spec-1), i,
                                 N_star_mono[i]))[i]
                    counter += 1
                if counter >= max_counter-2:
                    raise InputError(('Monoculture growth rate of species {i} '
                        'does not become negative with increasing N_{i}, '
                        'i.e. ``f_{i}(N_{i})``>0 for any N').format(i=i))
                N_star_mono[i] = brentq(lambda N: pars["f"](
                        np.insert(np.zeros(n_spec-1), i,N))[i],
                    0, N_star_mono[i])
            else:
                N_star_mono[i] = fsolve(pars["f"](
                        np.insert(np.zeros(n_spec-1), i,N))[i])
        # estimate that equilibrium density in each community is the 
        # monoculture equilibrium
        pars["N_star"] = N_star_mono*np.ones((n_spec, n_spec))    
    
    # c must be a positive real number
    if (np.any(~np.isfinite(pars["c"])) or np.any(pars["c"]<=0) 
            or pars["c"].dtype != float):
        warn("Some entries in pars['c'] were not positive real numbers."
             "These are replaced with 1")
        pars["c"] = np.real(pars["c"])
        pars["c"][pars["c"] <= 0] = 1
        pars["c"][~np.isfinite(pars["c"])] = 1
    
    for i in range(n_spec):
        # to set species i to 0
        ind = np.arange(n_spec) != i
        # solve for equilibrium, use equilibrium dens. of previous run
        N_pre,info,a ,b = fsolve(lambda N: pars["f"](np.insert(N,i,0))[ind],
                            pars["N_star"][i,ind], full_output = True,
                            xtol = xtol)
        
        # Check stability of equilibrium
        # Jacobian of system at equilibrium
        r = np.zeros((n_spec-1, n_spec-1))
        r[np.triu_indices(n_spec-1)] = info["r"].copy()
        jac = np.diag(N_pre).dot(info["fjac"].T).dot(r)
        # check whether we found equilibrium
        if np.amax(np.abs(info["fvec"]))>xtol:
            pars["equilibrium found with spec{} absent".format(i)] = N_pre
            pars["growth at found equilibrium"] = info["fvec"]
            try:
                pars["eigenvalues equilibrium"] = np.linalg.eigvals(jac)
            except np.linalg.LinAlgError:
                pass
            pars["fsolve output"] = info
            raise InputError("Not able to find resident equilibrium density, "
                        + "with species {} absent.".format(i)
                        + " Please provide manually via the `pars` argument")
        
        # check whether equilibrium is feasible, i.e. positive
        if not (np.all(N_pre>0) and np.all(np.isfinite(N_pre))):
            pars["equilibrium found with spec{} absent".format(i)] = N_pre
            pars["growth at found equilibrium"] = info["fvec"]
            try:
                pars["eigenvalues equilibrium"] = np.linalg.eigvals(jac)
            except np.linalg.LinAlgError:
                pass
            pars["fsolve output"] = info
            raise InputError("Found equilibrium is not feasible (i.e. N*>0), "
                        + "with species {} absent.".format(i)
                        + " Please provide manually via the `pars` argument")
        
        # check whether real part of eigenvalues is negative
        if max(np.real(np.linalg.eigvals(jac)))>0:
            pars["equilibrium found with spec{} absent".format(i)] = N_pre
            pars["growth at found equilibrium"] = info["fvec"]
            try:
                pars["eigenvalues equilibrium"] = np.linalg.eigvals(jac)
            except np.linalg.LinAlgError:
                pass
            pars["fsolve output"] = info
            raise InputError("Found equilibrium is not stable, "
                        + "with species {} absent.".format(i)
                        + " Please provide manually via the `pars` argument")        
            
        # save equilibrium density and invasion growth rate
        pars["N_star"][i] = np.insert(N_pre,i,0)
        pars["r_i"][i] = pars["f"](pars["N_star"][i])[i]
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
    # check for special cases first
    no_comp = np.isclose([NO_fun(pars,1, sp),
                          NO_fun(pars,1, sp[::-1])], [0,0])
    if no_comp.any():
        return special_case(no_comp, sp)
    
    sp = np.asarray(sp)
    
    def inter_fun(c):
        # equation to be solved
        NO_ij = np.abs(NO_fun(pars,c, sp))
        NO_ji = np.abs(NO_fun(pars,1/c,sp[::-1]))
        return NO_ij-NO_ji
    
    # use a generic numerical solver when `f` is not montone
    # potentially there are multiple solutions
    if not monotone_f: 
        print(pars["N_star"].shape)
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
    return c, 1/c # return c_i and c_j = 1/c_i

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
    
def NO_fun(pars,c, sp):
    # Compute NO for specis sp and conversion factor c
    f0 = pars["f"](switch_niche(pars["N_star"][sp[0]],sp))[sp[0]]
    fc = pars["f"](switch_niche(pars["N_star"][sp[0]],sp,c))[sp[0]]

    if f0 == fc:
        return np.sign(f0-pars["r_i"])[sp[0]]*np.inf
    return (f0-pars["r_i"][sp[0]])/(f0-fc)
    
def FD_fun(pars, c, sp):
    # compute the FD for species sp and conversion factor c
    f0 = pars["f"](switch_niche(pars["N_star"][sp[0]],sp))[sp[0]]
    fc = pars["f"](switch_niche(pars["N_star"][sp[0]],sp,c))[sp[0]]
    
    return fc/f0
    
def switch_niche(N,sp,c=0):
    # switch the niche of sp[1:] into niche of sp[0]
    N = N.copy()
    N[sp[0]] += np.nansum(c*N[sp[1:]])
    N[sp[1:]] = 0
    return N