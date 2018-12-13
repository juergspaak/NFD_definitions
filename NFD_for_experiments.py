"""
@author: J.W.Spaak
Numerically compute ND and FD for experimental data
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.optimize import brentq

try:
    from numerical_NFD import NFD_model
except ImportError:
    # in case this code is used in a submodule, import from the submodule
    from nfd_definitions.numerical_NFD import NFD_model

def NFD_experiment(N_star, time_exp1, dens_exp1, time_exp2, dens_exp2,
                     r_i, visualize = True):
    """Compute the ND and FD for two-species experimental data
    
    Compute the niche difference (ND), niche overlapp (NO), 
    fitnes difference(FD) and conversion factors (c). The 3 experiments to
    conduct are described in:
    The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
    DOI:
    
    Parameters
    -----------
    N_star: ndarray (shape = 2)
        Monoculture equilibrium density for both species
    time_exp1: array like
        Timepoints at which measurments of exp1 were taken in increasing order.
        If timepoints differ from species to species time_exp1 should be the
        union of all timepoints, and missing data should be indicated in
        exp1_data.
    dens_ex12: ndarray (shape = (2, len(time_exp1)))
        exp1_data[i,t] is the density of species i at time time_exp1[t] in exp1
        np.nan are allowed in case not all species have the same timeline for
        exp1
    time_exp2, dens_exp2: Equivalent to exp1
    r_i: ndarray (shape = 2)
        Invasion growth rate of both species
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
        only returned if visualize is True
    ax: Matplotlib axes
        only returned if visualize is True
        
    Literature:
    The unified Niche and Fitness definition, J.W.Spaak, F. deLaender
    DOI:
    """
    # combine all data into lists
    times = [None,time_exp1, time_exp2]
    exps = [None,dens_exp1, dens_exp2]
    
    # monoculture growth rate, assumes that growth rate is constant at
    # the beginning
    f0 = np.log(dens_exp1[:,1]/dens_exp1[:,0])/(time_exp1[1]-time_exp1[0])
    
    # per capita growth rate for both species in monoculture
    f = per_capita_growth(times, exps, N_star,f0)
    
    # compute the ND, FD etc. parameters
    pars = {"N_star": np.array([[0,N_star[1]],[N_star[0],0]]),
            "r_i": r_i}    
    pars = NFD_model(f, pars = pars, experimental = True)
    pars["f"] = f
    
    if visualize: # visualize results if necessary
        fig, ax = visualize_fun(f,times, exps, N_star, pars)
        return pars, fig, ax
    
    return pars
    
def per_capita_growth(times, exps, N_star,f0):
    """interpolate the per capita growth rate of the species
    
    times:
        Timepoints of measurments
    exps:
        Densities of the species at timepoints times
    N_star: float
        equilibrium density of the species
    f0: float
        monoculture growth rate
        
    Returns
    -------
    f: callable
        Per capita growth rate of species, fullfilling the differential
        equation dN/dt=N*f(N)
        Values below min(exps) are assumed to be f0, Values above max(exps)
        are assumed to be f(max(exps))
    """
    # percapita growth rates for each of the experiments separately
    dict_subf = {"f_exp{}_spec{}".format(i,j):
                dens_to_per_capita(times[i], exps[i][j])
                for i in [1,2] for j in [0,1]}
    
    # interpolation for datas between the two experiments   
    inter_data0 = np.array([
            [exps[1][0,-1], dict_subf["f_exp1_spec0"](exps[1][0,-1])],
            [N_star[0],0],
            [exps[2][0,-1], dict_subf["f_exp2_spec0"](exps[2][0,-1])]
            ])
    inter_data1 = np.array([
            [exps[1][1,-1], dict_subf["f_exp1_spec1"](exps[1][1,-1])],
            [N_star[1],0],
            [exps[2][1,-1], dict_subf["f_exp2_spec1"](exps[2][1,-1])]
            ])
    # quadrativ interpolation
    dict_subf["f_mid_spec0"] = ius(*inter_data0.T, k = 2)
    dict_subf["f_mid_spec1"] = ius(*inter_data1.T, k = 2)
    
    # per capita growth rate for each species, using different cases
    def f_spec(N,i):
        if N<exps[1][i,1]: # below minimum, use f0
            return f0[i]
        elif N<exps[1][i,-1]: # use values of exp1
            return dict_subf["f_exp1_spec"+str(i)](N)
        elif N<exps[2][i,-1]: # values between the two experiments
            return dict_subf["f_mid_spec" +str(i)](N)
        elif N<=exps[2][i,0]*0.99: # use values of exp2
            return dict_subf["f_exp2_spec"+str(i)](N)
        else: # above maximum
            return dict_subf["f_exp2_spec"+str(i)](exps[2][i,0]*0.99)
    
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
    
    return f
        
            
def dens_to_per_capita(time, dens, k = 3):
    # convert densities over time to per capita growth rate
    
    # remove nan's
    ind = np.isfinite(dens)
    time = time[ind]
    dens = dens[ind]
    
    # interpolate the data
    N_t = ius(time,dens,k=k)
    dNdt = N_t.derivative() # differentiate
    def per_capita(N):
        # search for t with N(t) = N
        try:
            t_N = brentq(lambda t: N_t(t)-N,time[0], time[-1])
        except ValueError:
            print(N,N_t(time[0]), N_t(time[1]))
            raise
        return dNdt(t_N)/N

    return per_capita

def visualize_fun(f,times, exps, N_star, pars):
    # visualize the fitted population density and the per capita growthrate
    fig, ax = plt.subplots(2,2, figsize = (9,9), sharey = "row", sharex ="row")
    
    # plot the densities over time
    # plot real data of exp1 for both species
    ax[0,0].scatter(times[1], exps[1][0], color = "black")
    ax[0,1].scatter(times[1], exps[1][1], color = "black")
    
    # plot real data of exp2 for both species
    ax[0,0].scatter(times[2], exps[2][0],facecolor = "none", color = "black")
    ax[0,1].scatter(times[2], exps[2][1],facecolor = "none", color = "black")
    
    time_1 = np.linspace(*times[1][[0,-1]], 100)
    time_2 = np.linspace(*times[2][[0,-1]], 100)
    
    # plot fitted data of exp1 for both species
    ax[0,0].plot(time_1, ius(times[1],exps[1][0])(time_1), label = "fit, exp1")
    ax[0,0].plot(time_2, ius(times[2],exps[2][0])(time_1), label = "fit, exp2")
    
    # plot fitted data of exp1 for both species
    ax[0,1].plot(time_1, ius(times[1],exps[1][1])(time_1), label = "fit, exp1")
    ax[0,1].plot(time_2, ius(times[2],exps[2][1])(time_1), label = "fit, exp2")
    
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
    
    
    # plot the fitted per capita growth rate
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