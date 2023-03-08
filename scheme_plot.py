import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

def do_all(r_i, **kwargs):
    """
    Given the invasion scheme plots the invsaion graph and finds end states
    
    Given all the feasible sub-communities for a given community model it
        
    computes the invasion graph and plots it.
    This function is solely a wrapper which calls the functions
    compute_graph and plot_scheme for convenience
    
    Parameters
    ----------
    r_i : numpy array, shape (n_equi, n_spec)
        Invasion scheme containing the invasion growth rates of the ``n_spec``
        species in the ``n_equi`` different equilibrium densities
        Non-stable equilibria may be included with np.nan rows
    **kwargs: Any keyword which can be passed to the function plot_scheme

    Returns
    -------
    
    invasion_graph: numpy array, shape (n_equi, n_equi), dtype bool
        Adjacency matrix of invasion graph
    cyclic: boolean
        Whether the invasion graph contains any cycles
    """
    
    # compute invasion graph
    inv_graph, cycles, permanent = compute_graph(r_i)
    # plot graph
    plot_scheme(r_i, **kwargs)
    
    return inv_graph, cycles, permanent 

def compute_graph(inv_scheme):
    """
    Computes the invasion graph for a given invasion scheme
    
    The invasion graph is a directed graph. A directed edge goes from vertix S
    to vertix T if and only if:
        For all i in S\T, i can not invade the state T
        For all j in T\S, j can invade in the state S
    Additionally it verifies whether the invasion graph contains any cycles

    Parameters
    ----------
    inv_scheme : numpy array, shape (n_equi, n_spec)
        Invasion scheme containing the invasion growth rates of the ``n_spec``
        species in the ``n_equi`` different equilibrium densities
        Only stable and feasible equilibria should be present

    Returns
    -------
    invasion_graph: numpy array, shape (n_equi, n_equi), dtype bool
        Adjacency matrix of invasion graph
    cycles: list
        List of all simple cycles in the invasion graph
    permanent: list
        List of the communities which are permanent
    """
    
    # create adjacency matrix of graph
    adj = np.full((2, len(inv_scheme), len(inv_scheme)), False, dtype = bool)
    # fill adjacency matrix accordint to the two rules
    for i in range(len(inv_scheme)):
        adj[0, i] = np.all((inv_scheme != 0) | (inv_scheme[i] >= 0), axis = -1)
        adj[1, i] = np.all((inv_scheme <= 0) | (inv_scheme[i] != 0), axis = -1)
    # only a directed link if both conditions are satisfied
    adj_real = np.all(adj, axis = 0)
    
    # remove self loops
    adj_real[np.arange(adj_real.shape[-1]),
             np.arange(adj_real.shape[-1])] = False
    
    A = adj_real.copy()
    # check for cycles, compute adj_real^len(adj_real)
    for i in range(2+int(np.log2(len(A)))):
        A = A.dot(A)
    
    if A.any():
        return adj_real, True, []
    
    # find all permanent subcommunities
    # remove non-present communities
    inv_scheme = inv_scheme[np.all(np.isfinite(inv_scheme), axis = 1)]
    # check for saturated communities
    
    saturated = np.where(np.all(inv_scheme<=0, axis = 1))[0]
    
    permanent = []
    for i in saturated:
        sp_present = inv_scheme[i] == 0
        
        # sub invasion scheme where only those species are present
        # remove equilibria where another species is present
        only_sp_present = np.all((inv_scheme == 0) <= sp_present, axis = 1)
        sub_invasion_scheme = inv_scheme[only_sp_present]
        # remove invasion growth rates of other species
        sub_invasion_scheme = sub_invasion_scheme[:,sp_present]
        
        # remove permanent subcommunity itself
        sub_invasion_scheme = sub_invasion_scheme[
                    ~np.all((sub_invasion_scheme == 0), axis = 1)]
        
        # check whether in each subcommunity at least one species can invade
        if np.all(np.any(sub_invasion_scheme>0, axis = 1), axis = 0):
            permanent.append(i)
    
    return adj_real, [], permanent

def plot_scheme(inv_scheme, axc = None, inv_graph = None,
                fs = 10, radius = 0.05, plot_non_stable = True,
                sp_names = None, plot_min_i_com = False):
    """
    PLots the invasion graph given the invasion scheme

    Parameters
    ----------
    inv_scheme : numpy array, shape (n_equi, n_spec)
        Invasion scheme containing the invasion growth rates of the ``n_spec``
        species in the ``n_equi`` different equilibrium densities
        Non-stable equilibria may be included with np.nan rows
    inv_graph: numpy array, shape (n_equi, n_equi), dtype bool, optional
        Adjacency matrix of invasion graph.
        Can be computed automatically if not given
    axc: axes into which the graph should be plotted
    fs: Fontsize of labeling for species
    plot_non_stable: boolean, default True
        Whether or not to also plot the non-existant subcommunities.
        Recommended to set to false for species rich communities (n>5)

    Returns
    -------
    None
    """    
    if axc is None:
        plt.figure()
        axc = plt.gca()
    
    if inv_graph is None:
        inv_graph, acyclic, id_permanent = compute_graph(inv_scheme)
    n = inv_scheme.shape[-1]
    
    # convert inv_scheme into indices
    id_inv_scheme = np.sum((inv_scheme == 0)*2**np.arange(n), axis = 1)
    id_saturated = id_inv_scheme[id_permanent]
    
    # plot all possible sub-communities or only stable ones?
    if plot_non_stable:
        indices = np.arange(2**n)
    else:
        indices = id_inv_scheme
    
    # locations on x and y axis
    x_locs, y_locs = np.empty((2,2**n))
    sp_present = np.empty(2**n, dtype = "object")
    richness = np.zeros(2**n)
    n_specs = np.arange(n+1)
    
    for i in n_specs:
        combs_sp = np.array(list(combinations(np.arange(n), i)))
        combs = np.sum(2**combs_sp, axis = 1).astype(int)
        stable_subcoms = combs[[(a in indices) for a in combs]]
        
        x_locs[stable_subcoms] = np.linspace(0,1, len(stable_subcoms) + 2)[1:-1]
        y_locs[stable_subcoms] = .2*np.linspace(-1,1,len(stable_subcoms))**2
        y_locs[stable_subcoms] -= np.mean(y_locs[stable_subcoms]) - i
        if i == 0:
            continue
        sp_present[combs] = [a for a in combs_sp]
        richness[combs] = i
    sp_present[0] = []
    
    # rescale y-location to fit into 0 and 1
    y_locs = y_locs/np.amax(y_locs)
    y_locs = (y_locs - 0.5)*0.9 + 0.5
    
    
    if sp_names is None:
        sp_names = n_specs
    elif (type(sp_names) == str) and (sp_names == "index_at_1"):
        sp_names = n_specs + 1
    
    # plot the vertices of the graph
    for i in indices:
        sp = sp_present[i]
        if len(sp) == 0:
            text = "-"
        else:
            order = int(np.sqrt(len(sp)))
            order += 1 if order**2 != len(sp) else 0
            text_array = np.empty((len(sp), 2), dtype = "object")
            text_array[:,0] = sp_names[sp]
            text_array[:,1] = ","
            text_array[order-1::order,1] = "\n"
            text_array[-1,-1] = ""
            text = "".join([str(i) for i in text_array.flatten()])
        if i in id_saturated:
            color = "purple"
        elif i in id_inv_scheme:
            color = "k"
        else:
            color = "red"
        axc.text(x_locs[i], y_locs[i], text
                 , ha = "center", va = "center",
                 fontsize = fs)
        circ = plt.Circle([x_locs[i], y_locs[i]], radius = radius,
                          edgecolor = color, facecolor = "white", transform = axc.axes.transAxes)
        axc.add_patch(circ)
        
    # should -i communities be highlighted?
    if plot_min_i_com:
        min_i_com = np.where(np.sum(inv_scheme>0, axis = 1)<=1)[0]
        for i in min_i_com:
            circ = plt.Circle([x_locs[i], y_locs[i]], radius = radius,
                          edgecolor = None, facecolor = "yellow",
                          transform = axc.axes.transAxes)
            axc.add_patch(circ)
    
    # plot the edges of the graph
    for k, i in enumerate(id_inv_scheme):
        for id_j, j in enumerate(inv_graph[k]):
            if j:    
                vector = np.array([[x_locs[id_inv_scheme[id_j]],
                                    y_locs[id_inv_scheme[id_j]]],
                                   [x_locs[i], y_locs[i]]])
                direction = vector[0]-vector[1]
                vector[0] -= direction/np.linalg.norm(direction)*radius
                diff_rich = richness[id_inv_scheme[id_j]] - richness[i]
                
                color = ["red", "orange", "blue"][int(np.sign(diff_rich)+1)]
                if diff_rich>1:
                    color = "lightgrey"
                axc.annotate("", xy = vector[0], xytext = vector[1],
                             arrowprops = {"arrowstyle": "->",
                                           "edgecolor": color,
                                           "alpha": 0.5},
                             zorder = 0)
    
    
    axc.set_yticks(np.linspace(min(y_locs), max(y_locs), int(max(richness[id_inv_scheme])+1)))
    axc.set_yticklabels(np.arange(int(max(richness[id_inv_scheme])+1)))

###############################################################################
# Lotka-Volterra related functions
     
def do_all_LV(A, mu = None, method = "old", **kwargs):
    """
    Computes the invasion scheme and plots graph for a LV community model
    
    Computes all feasible subcommunities and the corresponding invasion growth
    rates for a community model of the form
        1/N_i dN_i/dt = \mu_i - \sum_j a_ij N_j
        
    Next it computes the invasion graph and plots it.
    This function is solely a wrapper which calls the functions
    compute_scheme_for_LV, compute_graph and plot_scheme for convenience
    
    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates
    **kwargs: Any keyword which can be passed to the function plot_scheme

    Returns
    -------
    inv_scheme : numpy array, shape (n_equi, n_spec)
        Invasion scheme containing the invasion growth rates of the ``n_spec``
        species in the ``n_equi`` different equilibrium densities
        Non-stable equilibria may be included with np.nan rows
    invasion_graph: numpy array, shape (n_equi, n_equi), dtype bool
        Adjacency matrix of invasion graph
    cyclic: boolean
        Whether the invasion graph contains any cycles
    permanent: integer
        index of permanent end state of assembly
    pars : dictionary, keys are:
        ``ND``: numpy array, shape (n_spec,) contains the niche differences
        ``FD``: numpy array, shape (n_spec,) contains the fitness differences
        ``mu_i``: numpy array, shape (n_spec,) intrinsic growth rate
        ``r_i``: numpy array, shape (n_spec,), invasion growth rate
        ``eta_i``: numpy array, shape (n_spec,) no_niche growth rate
        ``c``: numpy array, shape (n_spec, n_spec), conversion factors
    """
    
    # compute invasion scheme
    r_i, equi_all = compute_scheme_for_LV(A, mu)
    # compute invasion graph
    inv_graph, cycles, permanent = compute_graph(r_i)
    # plot graph
    plot_scheme(r_i, **kwargs)
    
    pars = compute_NFD_LV(A, mu, method)
    
    return r_i, inv_graph, cycles, permanent, pars

def compute_scheme_for_LV(A, mu = None):
    """
    Computes the invasion scheme for a LV community model
    
    Computes all feasible subcommunities and the corresponding invasion growth
    rates for a community model of the form
        1/N_i dN_i/dt = \mu_i - \sum_j a_ij N_j
    

    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates

    Returns
    -------
    inv_scheme : numpy array, shape (n_equi, n_spec)
        Invasion scheme containing the invasion growth rates of the ``n_spec``
        species in the ``n_equi`` different equilibrium densities
        Non-stable equilibria may be included with np.nan rows
    """
    
    n = len(A)
    if mu is None:
        mu = np.ones(n)
    counter = 0
    equi_all = np.zeros(( 2**n, n))
    r_i = np.full(equi_all.shape, np.nan)
    # compute all possible equilibria
    for n_spec in range(1, n+1):
        combs = np.array(list(combinations(range(n), n_spec)))
        A_temp = np.empty((len(combs), n_spec, n_spec))
        mu_temp = np.empty((len(combs), n_spec))
        for i, comb in enumerate(combs):
            A_temp[i] = A[comb[:,np.newaxis], comb]
            mu_temp[i] = mu[comb]
        
        # solve equilibrium vectorized for speed
        equi_temp = np.linalg.solve(A_temp, mu_temp)
        
        for i, comb in enumerate(combs):
            counter += 1
            equi_all[counter, comb] = equi_temp[i]
    # community when all species are absent
    equi_all[0] = 0
    
    # compute invasion growth rates
    r_i = mu - np.sum(A*equi_all[:,np.newaxis], axis = -1)
    r_i = r_i[~np.any(equi_all<0, axis = 1)]
    
    # remove rounding errors
    r_i = np.round(r_i, 10)
    
    return r_i, equi_all

def equi_for_invasion(A, mu = None, remove_non_feasible = True,
                      ret_groups = False):
    """For each species computes the resident community at the end state
    
    Given the LV-community  given by
        1/N_i dN_i/dt = \mu_i - \sum_j a_ij N_j
    it computes the invasion graph and a possible end state of assembly.
    Given this end state of assembly it computes the resident community for all
    species at invasion.
    For the species not present at the end-state, the resident community is
    simply the end state of assembly.
    For the species present the resident community is the n-1 community of
    the end state of assembly. If this n-1 community does not exist then the
    resident community is np.nan.

    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates. Default is mu = 1
    remove_non_feasible : Bool, optional
        Whether non-feasible resident communities are replaced with np.nan.
        Recommended unless in debugging. The default is True.
    ret_groups : bool, optional
        Whether to return the different species groups along with the resident
        communities. If ''True'', then excl, w_n1 and wo_n1 are returned.
        The default is False.

    Raises
    ------
    RuntimeError
        If there are more than one end state to assembly

    Returns
    -------
    equi_invasion: numpy array, shape (n_spec, n_spec)
        equi_invaison[i] is the density of the resident community when species
        i invades
    excl: numpy array
        indices of species which are excluded from the end state of assembly.
        Only returned if ret_groups is True
    wo_n1: numpy array:
        Indices of species which are present at the end state of assembly, but
        do not have an n-1 community. For these species equi_invasion[i] is np.nan
        Only returned if ret_groups is True
    wo_n1: numpy array:
        Indices of species which are present at the end state of assembly, and
        do have an n-1 community. For these species equi_invasion[i] is np.nan
        Only returned if ret_groups is True

    """
    if mu is None:
        mu = np.ones(len(A))
    
    # compute equilibrium and invasion growth rates
    r_i, equi_all = compute_scheme_for_LV(A, mu)
    adj_real, cyclic, permanent = compute_graph(r_i)
    if len(permanent) != 1:
        raise RuntimeError("Feature not implmented yet")
    
    # species present
    pres = np.where(r_i[permanent[0]] == 0)[0]
    # species which are excluded
    excl = np.where(r_i[permanent[0]] != 0)[0]
    
    # compute the resident community for all invasions
    equi_invasion = np.zeros((len(mu), len(mu)))
    equi_invasion[excl[:,np.newaxis], pres] = np.linalg.solve(A[pres[:,np.newaxis], pres], mu[pres])
    
    # compute invasion densities for the other species
    for i in pres:
        sub_p = pres[pres != i]
        equi_invasion[i, sub_p] = np.linalg.solve(A[sub_p[:,np.newaxis], sub_p],
                                             mu[sub_p])
    
    if remove_non_feasible:
        # some equilibria might not be feasible
        equi_invasion[np.any(equi_invasion<0, axis = 1)] = np.nan
    
    if ret_groups:
        return equi_invasion, excl, pres[np.isnan(equi_invasion[pres,0])], pres[np.isfinite(equi_invasion[pres,0])]
    return equi_invasion

def __NFD_simple__(A, mu, equi_invasion):
    """Computes niche and fitness differences for an LV community model    
    
    This function should not be called directly, but rather is called internally
    by compute_NFD_LV

    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates. Default is mu = 1
    equi_invasion: numpy array, shape (n_spec, n_spec)
        equi_invaison[i] is the density of the resident community when species
        i invades. Optimally this is computed with ''equi_for_invasion""

    Returns
    -------
    pars : dictionary, keys are:
        ``ND``: numpy array, shape (n_spec,) contains the niche differences
        ``FD``: numpy array, shape (n_spec,) contains the fitness differences
        ``mu_i``: numpy array, shape (n_spec,) intrinsic growth rate
        ``r_i``: numpy array, shape (n_spec,), invasion growth rate
        ``eta_i``: numpy array, shape (n_spec,) no_niche growth rate
        ``c``: numpy array, shape (n_spec, n_spec), conversion factors
    """
    # compute conversion factor
    c_two = np.sqrt(np.abs(A/A.T/np.diag(A)[:,np.newaxis]*np.diag(A)))
    
    # invasion growth rate
    r_i = mu - np.diag(np.einsum("ij, nj->ni", A, equi_invasion))
    
    # no-niche growth rate
    eta_i = mu - np.diag(A)*np.diag(np.einsum("ij, nj->ni", c_two, equi_invasion))
    
    ND = (r_i - eta_i)/(mu - eta_i)
    F = -eta_i/(mu - eta_i)
    
    pars = {"ND": ND, "F": F,
            "mu_i": mu, "r_i": r_i, "eta_i": eta_i, "c": c_two}
    
    return pars

def change_focus_LV(A, mu, focus, non_foc):
    """
    Change focus of LV community model
    
    For a given LV community model and a set of focal species, this computes
    a new LV community model where the non focal species are treated
    similar to resaources and are assumed to be at equilibrium.

    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates. Default is mu = 1
    focus : numpy array, dtype int
        Indices of focal species
    non_foc :numpy array, dtype int
        INdices of non-focal species

    Returns
    -------
    A_new: numpy array, shape (len(focus), len(focus))
        Community interaction matrix for focal species
    mu_new: numpy array, shape (n_spec), optional
        Intrinsic growth rates for focal species

    """
    A_new = (A[focus[:,np.newaxis], focus] -
            A[focus[:,np.newaxis], non_foc].dot(
                np.linalg.inv(A[non_foc[:,np.newaxis], non_foc])
                ).dot(A[non_foc[:,np.newaxis], focus]))
    mu_new = (mu[focus] - 
              A[focus[:,np.newaxis], non_foc].dot(
                  np.linalg.inv(A[non_foc[:,np.newaxis], non_foc])
                  .dot(mu[non_foc])))
    return A_new, mu_new
    

def compute_NFD_LV(A, mu = None, method = "old"):
    """Computes niche and fitness differences for an LV community model    

    Parameters
    ----------
    A: numpy array, shape (n_spec, n_spec)
        Community interaction matrix
    mu: numpy array, shape (n_spec), optional
        Intrinsic growth rates. Default is mu = 1
    method: "old" or "new", default = "old"
        Affects how niche and fitness differences for present species should
        be computed, specifically the conversion factor ``c_ij``. IF old, then
        the conversion factors are computed based on the interaction matix A.
        If "new", then the conversion factors are computed based on a new 
        interaction matrix A_new, which is obtained by treating the species
        without a n-1 community similar to resources.

    Returns
    -------
    pars : dictionary, keys are:
        ``ND``: numpy array, shape (n_spec,) contains the niche differences
        ``FD``: numpy array, shape (n_spec,) contains the fitness differences
        ``mu_i``: numpy array, shape (n_spec,) intrinsic growth rate
        ``r_i``: numpy array, shape (n_spec,), invasion growth rate
        ``eta_i``: numpy array, shape (n_spec,) no_niche growth rate
        ``c``: numpy array, shape (n_spec, n_spec), conversion factors
    """
    if mu is None:
        mu = np.ones(len(A))
    # compute the invasion growh
    equi_invasion, excl, focus, non_foc = equi_for_invasion(A, mu, ret_groups=True)
    
    pars = __NFD_simple__(A, mu, equi_invasion)
    if method == "old":
        return pars

    # change focus of species
    A_foc, mu_foc = change_focus_LV(A, mu)
    pars_new = __NFD_simple__(A_foc, mu_foc)
    
    for key in pars.keys():
        if key == "c":
            continue
        else:
            pars[key][focus] = pars_new[key]
    pars["c_new"] = pars_new["c"]
    return pars


if __name__ == "__main__":
    ###########################################################################
    # test correctness of computation for simple niche and fitness differences
    n_spec = 3
    A = np.random.normal(0.2,0.1, (n_spec, n_spec))
    np.fill_diagonal(A, np.random.normal(1, 0.05, n_spec))
    mu = np.ones(len(A))
    equi_invasion = equi_for_invasion(A, mu)
    pars_new = compute_NFD_LV(A)

    from nfd_definitions.numerical_NFD import NFD_model
    pars = NFD_model(lambda N: mu - A.dot(N), n_spec = 3)
    
    for key in ["ND", "F", "c"]:
        if not np.allclose(pars[key], pars_new[key]):
            raise RuntimeError("Automatic and manual computations do not match for {}".format(key))
            
    ###########################################################################
    # test correctness of change of focus
    # equilibrium densities should not change
    
    n_spec = 5
    A = np.random.normal(0.2,0.1, (n_spec, n_spec))
    np.fill_diagonal(A, 1)
    mu = np.ones(len(A))
    focus = np.arange(3)
    non_foc = np.array([3,4])
    
    A_foc, mu_foc = change_focus_LV(A, mu, focus, non_foc)
    if not np.allclose(np.linalg.solve(A, mu)[focus],
                       np.linalg.solve(A_foc, mu_foc)):
        raise
