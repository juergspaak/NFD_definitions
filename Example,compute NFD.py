"""
@author: J.W.Spaak
Example how to compute the ND and FD for a given differential equation setting
"""

import numpy as np

try:
    from numerical_NFD import NFD_model
except ImportError:
    # in case this code is used in a submodule, import from the submodule
    from nfd_definitions.numerical_NFD import NFD_model

# create the differential equation system
n_spec = 2 # number of species in the system
np.random.seed(6) # set random seed for reproduce ability

# Lotka-Volterra model
A = np.random.uniform(0,1,(n_spec,n_spec)) # interaction matrix
np.fill_diagonal(A,np.random.uniform(1,2,n_spec)) # to ensure coexistence
mu = np.random.uniform(1,2,n_spec) # intrinsic growth rate
def test_f(N):
    return mu - np.dot(A,N)

# compute relevant parameters with software
pars = NFD_model(test_f, n_spec)
ND, NO, FD, c = pars["ND"], pars["NO"], pars["FD"], pars["c"]
# manualy check results for the two species case
# see appendix for proof of correctness
NO_check = np.sqrt(np.array([A[0,1]*A[1,0]/A[1,1]/A[0,0],
                       A[0,1]*A[1,0]/A[1,1]/A[0,0]]))
ND_check = 1-NO
FD_check = 1- mu[::-1]/mu*np.sqrt(np.array([A[0,1]*A[0,0]/A[1,0]/A[1,1],
                       A[1,0]*A[1,1]/A[0,1]/A[0,0]]))
c_check = np.sqrt(np.array([A[0,1]*A[1,1]/A[1,0]/A[0,0],
                      A[1,0]*A[0,0]/A[0,1]/A[1,1]]))

# precision of output
prec = 4
print("Results of two species case:\n")
print("\t Software\t\t Manual check\t\t Rel. Difference\n")
print("ND:\t", np.round(ND,prec), "\t", np.round(ND_check,prec),
      "\t", np.abs(ND-ND_check)/ND_check)
print("NO:\t", np.round(NO,prec), "\t", np.round(NO_check,prec),
      "\t", np.abs(NO-NO_check)/NO_check)
print("FD:\t", np.round(FD,prec), "\t", np.round(FD_check,prec),
      "\t", np.abs(FD-FD_check)/FD_check)
print("c:\t", np.round(c[[0,1],[1,0]],prec), "\t", np.round(c_check,prec),
      "\t", np.abs(c[[0,1],[1,0]]-c_check)/c_check)


###############################################################################
# Switching to multispecies case
# create the differential equation system
n_spec = 10 # nuber of species in the system

# Lotka-Volterra model
A = np.random.uniform(0,1,(n_spec,n_spec)) # interaction matrix

# to ensure coexistence increase diagonal values
np.fill_diagonal(A,np.random.uniform(n_spec,n_spec+1,n_spec)) 
mu = np.random.uniform(1,2,n_spec) # intrinsic growth rate
def test_f(N, mu, A):
    return mu - np.dot(A,N)

# how to pass additional arguments to NFD_model
pars = NFD_model(test_f, n_spec, args = (mu, A))
ND_m, NO_m, FD_m, c_m = pars["ND"], pars["NO"], pars["FD"], pars["c"]

NO_check_m = np.empty(n_spec)
FD_check_m = np.empty(n_spec)
for i in range(n_spec):
    denominator = 0
    numerator = 0
    for j in range(n_spec):
        if i==j:
            continue
        numerator += pars["N_star"][i,j]*A[i,j]
        denominator += pars["N_star"][i,j]*np.sqrt(A[i,j]/A[j,i]*A[i,i]*A[j,j])
    NO_check_m[i] = numerator/denominator
    FD_check_m[i] = 1-denominator/mu[i]
    
# printing layout is optimized for 6 species
def print_function(var, var_check, name):
    rel_diff = np.round(np.abs(var-var_check)/var_check,prec)
    var = np.round(var,prec)
    var_check = np.round(var_check,prec)
    print(name+":\t", var[:2],"\t", var_check[:2], "\t", rel_diff[:2])
    for i in range(1,len(var)//2):
        ind = [2*i,2*i+1]
        print("\t", var[ind],"\t", var_check[ind], "\t", rel_diff[ind])
    if len(var)%2==1:
        print("\t", var[-1],"\t\t", var_check[-1], "\t\t", rel_diff[-1])
    print()
        

print("\n\nResults of multi species case:\n")
print("\t Software\t\t Manual check\t\t Rel. Difference\n")
print_function(NO_m, NO_check_m, "NO")
print_function(FD_m, FD_check_m, "FD")
