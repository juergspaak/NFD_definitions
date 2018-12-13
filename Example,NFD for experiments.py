"""
@author: J.W.Spaak
Example how to compute the ND and FD for a given differential equation setting
"""

import numpy as np

from scipy.integrate import odeint
import matplotlib.pyplot as plt
try:
    from numerical_NFD import NFD_model
    from NFD_for_experiments import NFD_experiment
except ImportError:
    # in case this code is used in a submodule, import from the submodule
    from NFD_definitions.numerical_NFD import NFD_model
    from NFD_definitions.NFD_for_experiments import NFD_experiment

###############################################################################
# First create experimental dataset
###############################################################################
# create the differential equation system
n_spec = 2 # number of species in the system

# Lotka-Volterra model
A = np.random.uniform(0.5,1,(n_spec,n_spec)) # interaction matrix
np.fill_diagonal(A,np.random.uniform(1,2,n_spec)) # to ensure coexistence
mu = np.random.uniform(1,2,n_spec) # intrinsic growth rate
def test_f(N,t = 0):
    return mu - np.dot(A,N)

# create dataset for monoculture
time_exp1 = time_exp2 = np.linspace(0,10,20)
dens_exp1, dens_exp2 = np.empty((2,2,len(time_exp1)))

dens_exp1[0] = odeint(lambda N,t: N*test_f(N),[1e-3,0],time_exp1)[:,0]
dens_exp1[1] = odeint(lambda N,t: N*test_f(N),[0,1e-3],time_exp1)[:,1]

dens_exp2[0] = odeint(lambda N,t: N*test_f(N),[2*dens_exp1[0,-1],0]
                            ,time_exp1)[:,0]
dens_exp2[1] = odeint(lambda N,t: N*test_f(N),[0,2*dens_exp1[1,-1]]
                            ,time_exp1)[:,1]
N_star = mu/np.diag(A)
r_i = mu-N_star[::-1]*A[[0,1],[1,0]]


###############################################################################
# compute NFD given experimental data
# set visualize to True to have visual confirmation that interpolation and
# differentiation are reasonable
pars, fig, ax = NFD_experiment(N_star, time_exp1, dens_exp1,
                    time_exp2, dens_exp2, r_i,  visualize = True)
plt.show()

print("\n\n")

###############################################################################
# check results with computation for model
pars_model = NFD_model(test_f)
prec = 4
ND, FD, c = pars["ND"], pars["FD"], pars["c"]
ND_check, FD_check, c_check = pars_model["ND"],  pars_model["FD"], pars_model["c"]
print("Compare resutls from model and experimental:\n")
print("\t Experimental\t\t Model check\t\t Rel. Difference\n")
print("ND:\t", np.round(ND,prec), "\t", np.round(ND_check,prec),
      "\t", np.abs(ND-ND_check)/ND_check)
print("FD:\t", np.round(FD,prec), "\t", np.round(FD_check,prec),
      "\t", np.abs(FD-FD_check)/FD_check)
print("c:\t", np.round(c[[0,1],[1,0]],prec), "\t", 
      np.round(c_check[[0,1],[1,0]],prec),
      "\t", np.abs(c[[0,1],[1,0]]-c_check[[0,1],[1,0]])/c_check[[0,1],[1,0]])