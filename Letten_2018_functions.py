import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.integrate import solve_ivp

data = pd.read_csv("Letten_2018.csv")

D = 0.033521
# check data for correctness
R_star = np.round(D*data["K"]/(data["mu"] - D), 4)

low = 10
med = 30
high = 50
# maximum growth rate
mu = np.array([data.loc[data.sucrose == low, "mu"].values,
              data.loc[data.sucrose == high, "mu"].values,
              data.loc[data.sucrose == med, "mu"].values])
# halfsaturation constant
K = np.array([data.loc[data.sucrose == low, "K"].values,
              data.loc[data.sucrose == high, "K"].values,
              data.loc[data.sucrose == med, "K"].values])
# nutrient quota
Q = data.loc[data.sucrose == 10, "Q"].values
# dilution per migration
d = 0.2
# resource at new location
S = 2.683
# number of species
n_spec = 4

def growth(t, N, case):
    X = N[:4]
    R = max(N[-1],0)
    
    dX_dt = mu[case]*R*X/(K[case] + R)
    dR_dt = -np.sum(Q*dX_dt)
    return np.append(dX_dt, dR_dt)
    
def dilution(dens):
    R = max(dens[-1],0)
    return d*dens[:4], (1-d)*(S-R)

def simulate_densities(X_start, n_cycles, n_pres = 1, R_start = S):
    R_start = R_start # initial concentration of resources is maximal
    time = np.arange(0,2*24+1e-5, n_pres) # two days, every hour
    densities = np.empty((n_spec+1, n_cycles, 2, len(time)))
    time_all = np.empty((n_cycles, 2, len(time)))
    for i in range(n_cycles):
        for case in range(2):
            sol = solve_ivp(growth, time[[0,-1]], np.append(X_start, R_start),
                            t_eval = time, args = (case,))
            densities[:,i,case] = sol.y
            X_start, R_start = dilution(sol.y[:,-1])
            time_all[i, case] = time + (i*2 + case)*48
    return densities.reshape((n_spec +1, -1)), time_all.reshape(-1)

def get_equilibrium_densities(pres):
    X_start = np.zeros(n_spec)
    X_start[pres] = 100
    
    densities, time_all = simulate_densities(X_start, 100, n_pres = 2)
    X_start, R_start = dilution(densities[:,-1])
    # compute high precicsion
    dens_precise, time_precise = simulate_densities(X_start, 10, n_pres = 0.1,
                                                    R_start = R_start)
    
    return dens_precise, time_precise


def decomposition(dens, time):
    case_const = np.full(time.shape, 2) # mean between conditions
    # create randomized conditions
    case_rand = np.random.uniform(0,1, time.shape) 
    case_rand[case_rand<np.median(case_rand)] = 0
    case_rand[case_rand!=0] = 1
    case_rand = case_rand.astype("int")
    
    case_real = np.empty(time.shape)
    case_real = case_real.reshape((-1, int(time[-1])//48))
    case_real[:,::2] = 0
    case_real[:,1::2] = 1
    case_real = case_real.flatten(order = "F")
    case_real = case_real.astype("int")
    
    dens_real = dens.copy()
    dens_real[dens_real == 0] = 1e-3
    dens_const = dens_real.copy()
    dens_const[-1] = np.mean(dens_const[-1])
    dens_rand = dens_real.copy()
    dens_rand[-1] = dens_rand[-1,
                    np.argsort(np.random.uniform(0,1,dens_rand[-1].shape))]
    
    r_is = {}
    r_is["const"] = np.array([growth(0, dens_const[:,i], case_const[i])/dens_const[:,i]
                          for i in range(len(time))]).T
    
    r_is["case"] = np.array([growth(0, dens_const[:,i], case_real[i])/dens_const[:,i]
                          for i in range(len(time))]).T
    
    r_is["res"] = np.array([growth(0, dens_real[:,i], case_const[i])/dens_real[:,i]
                          for i in range(len(time))]).T
    
    r_is["rand"] = np.array([growth(0, dens_rand[:,i], case_rand[i])/dens_rand[:,i]
                          for i in range(len(time))]).T
    
    r_is["real"] = np.array([growth(0, dens_real[:,i], case_real[i])/dens_real[:,i]
                          for i in range(len(time))]).T
    
    for key in r_is.keys():
        r_is[key] = np.mean(r_is[key], axis = 1)[:n_spec] + np.log(d)/48
        
    return r_is


if __name__ == "__main__":
    
    dens, time = get_equilibrium_densities([2,3])
    plt.plot(time, dens.T)

    ingr = decomposition(dens, time)
    print(ingr["real"])
