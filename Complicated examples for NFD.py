"""
@author: J.W.Spaak
In this example we show all possible pitfalls while using the NFD computation

The example serves as a mathematical illustration, the model is not
necessarily biologically meaningful. 
"""

import numpy as np
import matplotlib.pyplot as plt
try:
    from numerical_NFD import NFD_model
except ImportError:
    # in case this code is used in a submodule, import from the submodule
    from NFD_definitions.numerical_NFD import NFD_model

# to run without error change this variable to True
no_error = True

###############################################################################
# we first define the percapita growth rates
# the (monoculture) percapita growth rates are polynomials of degree 3
# figure 1 shows the per capita growth rate for the species in monoculutre
x = np.array([[0,3,6,9.0],
              [0,1,2,3]])
y = np.array([[2,-0.1,1,0],
              [2,1,1.5,0]])

coef = np.array([np.polyfit(x[0],y[0],3), np.polyfit(x[1],y[1],3)]).T

# interspecific competition is linear and chosen s.t. invasion growth rate 1.5
roots = np.array([np.real(np.roots(coef[:,0])[-1]),x[1,3]])
r_i = np.array([1.5,1.65])
comp = (coef[-1]-r_i)/roots[::-1]

# the percapita growth rate of the species
# the argument `coef` for demonstration of the args argument
def full_example(N, coef = coef):
    return coef[0]*N**3 +coef[1]*N**2 +coef[2]*N +coef[3] -comp*N[...,::-1]

# plot the funciton in monoculture
steps = 1000
end = 10
N_mono = np.zeros((steps,2))
N_mono[:,0] = np.linspace(0,end,steps)

percapita_mono = np.array([full_example(N_mono, coef)[:,0], 
                           full_example(N_mono[:,::-1], coef)[:,1]])

# plot the per capita growth rate in monoculture
plt.figure()
plt.plot(N_mono[:,0],percapita_mono.T, label = "species")
plt.grid()
plt.axis([0,end,-2,3])
plt.legend(["species 1", "species 2"])
plt.xlabel("Density (N)", fontsize = 14)
plt.ylabel(r"per capita growth rate $f_i(N_i,0)$", fontsize = 14)
plt.title("Per capita growth rate in monoculture")
plt.show()
print("Note that species 1 has two stable monoculture equilibria.")
print("Species 2 has a local minima (not equilibrium) at 1.",
      "The algorithm can't find equilibrium density in this case.",
      "\nWe have to provide a beter starting estimate for",
      "the equilibrium density (default is 1).")
print("Finally the f are not monotone, which has to be specified by:")
print("monotone_f = False\n\n\n")


if not no_error:
    print('\033[31m'+
        "Error created for illustration, set no_error to True\n\n\n")
    # will result in InputError, as fsolve can't find equilibrium density
    NFD_model(full_example, args= (coef,))
       
# Pass estimates of the equilibrium densities via the pars argument
# These only have to be estimates, not the actual equilibrium values

# Note: N_star_approx[0] is the equilibrium density of species 1, as species 0
# is absent
N_star_approx = np.array([[0,roots[1]],
                          [2,0]], dtype = float)
# compute ND etc.
pars = NFD_model(full_example, pars = {"N_star": N_star_approx.copy()},
                monotone_f = False, args = (coef,))


###############################################################################
# show potential other solutions to ND1 = ND2

c = 5**np.linspace(-3,1,101)
Nc = np.zeros((len(c),2))
Nc[:,0] = c*pars["N_star"][0,1]
NO_0 = (pars["f0"][0]-pars["r_i"][0])/(pars["f0"][0]-full_example(Nc)[:,0])


Nc = np.zeros((len(c),2))
Nc[:,1] = 1/c*pars["N_star"][1,0]
NO_1 = (pars["f0"][1]-pars["r_i"][1])/(pars["f0"][1]-full_example(Nc)[:,1])

# plot ND_1 = ND_2
plt.figure()
plt.title(r"Solve $ND_1 = ND_2$")
plt.xlabel(r"$c_1$", fontsize = 14)
plt.ylabel(r"$ND_i(c_1)$", fontsize = 14)

plt.plot(c,1-NO_0, label = r"$ND_1$")
plt.plot(c,1-NO_1, label = r"$ND_2$")
plt.grid()
plt.axis([0,4,-0.5,1.5])
plt.legend(fontsize = 14)
plt.show()

print("To find ND we have to solve ND1 = ND2.")
print("This equation has several solutions (as f is not monoton).")
print("Default will assume c=1 and solve the equation.")
print("To find other solutions we can pass the key 'c' to pars.\n\n\n")

###############################################################################
# Now we want to find all possible solutions

# compute ND etc.
c = np.ones((2,2))
c[0,1] = 0.5
pars1 = NFD_model(full_example, monotone_f = False,
                 pars = {"N_star": N_star_approx.copy(),"c":c})
c[0,1] = 3
pars2 = NFD_model(full_example, monotone_f = False,
                 pars = {"N_star": N_star_approx.copy(),"c":c})

if not no_error:
    print('\033[31m'+
        "Error created for illustration, set no_error to True\n\n\n")
    N_star_approx[1,0] = 3
    # will result in InputError, as equilibrium is not stable
    NFD_model(full_example, monotone_f = False,
             pars = {"N_star": N_star_approx.copy(),"c":c})

N_star_approx[1,0] = 9
# in this case we have to pass a starting value of c as an estimate aswell
# furthermore there's only one c
pars_new_eq = NFD_model(full_example, monotone_f = False,
                 pars = {"N_star": N_star_approx.copy(),"c":c})

plt.figure()
x = np.linspace(1e-5,1-1e-5,100)
plt.plot(x,x/(1-x), "red", label = "Coex. boundary")
plt.axis([0,1,-1,2])

plt.plot([pars["ND"][0],pars1["ND"][0], pars2["ND"][0]],
         [-pars["FD"][0],-pars1["FD"][0], -pars2["FD"][0]],'bo',
         label = "Species 1")
plt.plot(x,(x-r_i[0]/coef[-1,0])/(1-x), ":b")

plt.plot([pars["ND"][1],pars1["ND"][1], pars2["ND"][1]],
         [-pars["FD"][1],-pars1["FD"][1], -pars2["FD"][1]],'go',
         label = "Species 2")
plt.plot(x,(x-r_i[1]/coef[-1,1])/(1-x), ":g")

plt.plot(pars_new_eq["ND"][0], -pars_new_eq["FD"][0], 'b^')
plt.plot(pars_new_eq["ND"][1], -pars_new_eq["FD"][1], 'g^')

plt.legend(loc = "upper left")
plt.title("Possible ND and FD choices")
plt.xlabel("ND")
plt.ylabel(r"$-FD$")
plt.show()

print("This system has 4 different (correct) values for ND and FD per species.")
print("3 (dots) coming from the monoculture equilibria N1 ~ 2, N2 ~ 3.")
print("The 4. (triangle) comes from the different monoculture equilibrium N1 ~ 9.")
print("The ND-FD decomposition steming from the same invasion growth rate,"
      + " i.e. from the same resident equilibria, ")
print("will lie on the line -FD = ND/(1-ND)-r_i/f_i(0,0)/(1-ND)")