"""
@author: J.W.Spaak
Numerically compute ND and FD for a model
"""

import numpy as np
from scipy.optimize import brentq, fsolve
from warnings import warn

def NFD_model(f, n_spec = 2, args = (), monotone_f = True, pars = None,
             experimental = False, from_R = False):
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
        Any extra arguments to `f
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
    force: boolean, default False
        If True the parameters given in pars are assumed to be exact and are
        not checked. Should be set to True when used in combination with
        experimental data.
    from_R: boolean, default False
        Set to True if function is called via R by reticulate package.
        Converts types of f and equilibria.
        
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
        
    Literature:
    The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
    DOI: 10.1101/482703
    """ 
    if from_R:
        if n_spec-int(n_spec) == 0:
            n_spec = int(n_spec)
        else:
            raise InputError("Number of species (`n_spec`) must be an integer")
        fold = f

        def f(N, *args):
            # translate dataframes, matrices etc to np.array
            return np.array(fold(N,*args)).reshape(-1)  
    # check input on correctness
    monotone_f = __input_check__(n_spec, f, args, monotone_f)
    
    if experimental:
        if not ("c" in pars.keys()):
            pars["c"] = np.ones((n_spec, n_spec))
        if not ("r_i" in pars.keys()):
            pars["r_i"] = np.array([f(pars["N_star"][i], *args)[i] 
                        for i in range(n_spec)])
        pars["f"] = lambda N: f(N, *args)
    else:
        # obtain equilibria densities and invasion growth rates    
        pars = preconditioner(f, args,n_spec, pars)                  
    # list of all species
    l_spec = list(range(n_spec))
    # compute conversion factors
    c = np.ones((n_spec,n_spec))
    for i in l_spec:
        for j in l_spec:
            if i>=j: # c is assumed to be symmetric, c[i,i] = 1
                continue
            c[i,j] = solve_c(pars,[i,j], monotone_f[i] and monotone_f[j])
            c[j,i] = 1/c[i,j]

    # compute NO and FD
    NO = np.empty(n_spec)
    FD = np.empty(n_spec)
    
    for i in l_spec:
        # creat a list with i at the beginning [i,0,1,...,i-1,i+1,...,n_spec-1]
        sp = np.array([i]+l_spec[:i]+l_spec[i+1:])
        # compute NO and FD
        NO[i] = NO_fun(pars, c[i, sp[1:]], sp)
        FD[i] = FD_fun(pars, c[i, sp[1:]], sp)
    
    # prepare returning values
    pars["NO"] = NO
    pars["ND"] = 1-NO
    pars["FD"] = FD
    pars["c"] = c
    return pars
  
def __input_check__(n_spec, f, args, monotone_f):
    # check input on (semantical) correctness
    if not isinstance(n_spec, int):
        raise InputError("Number of species (`n_spec`) must be an integer")
    
    # check whether `f` is a function and all species survive in monoculutre
    try:
        f0 = f(np.zeros(n_spec), *args)
        if f0.shape != (n_spec,):
            raise InputError("`f` must return an array of length `n_spec`")   
    except TypeError:
        print("function call of `f` did not work properly")
        raise
    except AttributeError:
        fold = f
        f = lambda N, *args: np.array(fold(N, *args))
        f0 = f(np.zeros(n_spec), *args)
        warn("`f` does not return a proper `np.ndarray`")
        
    if min(f0)<=0 or (not np.all(np.isfinite(f0))):
        raise InputError("All species must have positive monoculture growth"
                    +"i.e. `f(0)>0`. Especially this value must be defined")
    
    # broadcast monotone_f if necessary
    return np.logical_and(monotone_f, np.full(n_spec, True, bool))
        
class InputError(Exception):
    pass
        
def preconditioner(f, args, n_spec, pars):
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
    pars_def = {"N_star": np.ones((n_spec, n_spec), dtype = "float"),
                "c": np.ones((n_spec,n_spec)),
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
    pars["f"] = lambda N: f(N,*args)
    
    for i in range(n_spec):
        # to set species i to 0
        ind = np.arange(n_spec) != i
        # solve for equilibrium, use equilibrium dens. of previous run
        N_pre,info,a ,b = fsolve(lambda N: pars["f"](np.insert(N,i,0))[ind],
                            pars["N_star"][i,ind], full_output = True)
        
        # check whether we found equilibrium
        if np.amax(np.abs(info["fvec"])/N_pre)>1e-10:
            raise InputError("Not able to find resident equilibrium density, "
                        + "with species {} absent.".format(i)
                        + " Please provide manually via the `pars` argument")
        
        # Check stability of equilibrium
        # Jacobian of system at equilibrium
        r = np.zeros((n_spec-1, n_spec-1))
        r[np.triu_indices(n_spec-1)] = info["r"]
        jac = -np.diag(N_pre).dot(info["fjac"]).dot(r)

        # check whether real part of eigenvalues is negative
        if max(np.real(np.linalg.eigvals(jac)))<0:
            raise InputError("Found equilibrium is not stable, "
                        + "with species {} absent.".format(i)
                        + " Please provide manually via the `pars` argument")
            
        # check whether equilibrium is feasible, i.e. positive
        if not (np.all(N_pre>0) and np.all(np.isfinite(N_pre))):
            raise InputError("Found equilibrium is not feasible (i.e. N*>0), "
                        + "with species {} absent.".format(i)
                        + " Please provide manually via the `pars` argument")
        
            
        # save equilibrium density and invasion growth rate
        pars["N_star"][i] = np.insert(N_pre,i,0)
        pars["r_i"][i] = pars["f"](pars["N_star"][i])[i]
        pars["f0"] = pars["f"](np.zeros(n_spec))
    return pars
    
def solve_c(pars, sp = [0,1], monotone_f = True):
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
    sp = np.asarray(sp)
    
    def inter_fun(c):
        # equation to be solved
        NO_ij = np.abs(NO_fun(pars,c, sp))
        NO_ji = np.abs(NO_fun(pars,1/c,sp[::-1]))
        return NO_ij-NO_ji
    
    # use a generic numerical solver when `f` is not montone
    # potentially there are multiple solutions
    if not monotone_f: 
        c = fsolve(inter_fun,pars["c"][sp[0],sp[1]])[0]
        if np.abs(inter_fun(c))>1e-10:
            raise ValueError("Not able to find c_{}^{}.".format(*sp) +
                "Please pass a better guess for c_i^j via the `pars` argument")
        return c
        
    # if `f` is monotone then the solution is unique, findint it with a more
    # robust method
        
    # find interval for brentq method
    a = pars["c"][sp[0],sp[1]]
    # find which species has higher NO for c0
    direction = np.sign(inter_fun(a))
    
    if direction == 0: # starting guess for c is correct
        return a
    fac = 2**direction
    if not np.isfinite(direction):
        raise InputError("function `f` seems to be returning nonfinite values")
    b = a*fac
    # change searching range to find c with changed size of NO
    while np.sign(inter_fun(b)) == direction:
        a = b
        b *= fac
    
    # solve equation
    try:
        return brentq(inter_fun,a,b)
    except ValueError:
        raise ValueError("f does not seem to be monotone. Please run with"
                         +"`monotone_f = False`")
    
def NO_fun(pars,c, sp):
    # Compute NO for specis sp and conversion factor c
    f0 = pars["f"](switch_niche(pars["N_star"][sp[0]],sp))[sp[0]]
    fc = pars["f"](switch_niche(pars["N_star"][sp[0]],sp,c))[sp[0]]
    
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