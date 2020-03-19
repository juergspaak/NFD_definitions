

# @author: J.W.Spaak
# Example how to compute the ND and FD for a given differential equation setting

# The script is originally written in python, we run python code from within R
# make sure that python is installed on the system used
# This is not part of the install.packages("reticulate") command!
library(reticulate)
# set your working directory
setwd("C:/Users/jspaak/Dropbox/Doktorat/Projects/P5_Fitness and niche differences/3_Programs")

# loads the relevant python code
source_python("numerical_NFD.py")

# create the differential equation system
n_spec <- 2 # number of species in the system, must be an integer
set.seed(0) # set random seed for reproduce ability

# Lotka-Volterra model
A <- matrix(runif(n_spec^2,0,1), n_spec, n_spec) # interaction matrix
diag(A) <- runif(n_spec, 1,2) # to ensure coexistence
mu <- runif(n_spec,1,2) # intrinsic growth rate
test_f <- function(N){
  return(mu - A%*%N)
}

# compute relevant parameters with python
# the parameter `from_R = TRUE` changes data types from R to python
pars <- NFD_model(test_f, n_spec, from_R = TRUE)
ND <- pars$ND
NO <- pars$NO
FD <- pars$FD
c <- pars$c

# manualy check results for the two species case
# see appendix for proof of correctness
NO_check = sqrt(A[1,2]*A[2,1]/(A[1,1]*A[2,2]))*c(1,1)
ND_check = 1-NO_check
FD_check = 1- rev(mu)/mu*sqrt(c(A[1,2]*A[1,1]/A[2,1]/A[2,2],
                                A[2,1]*A[2,2]/A[1,2]/A[1,1]))
c_check = sqrt(c(A[1,2]*A[2,2]/A[2,1]/A[1,1],
                 A[2,1]*A[1,1]/A[1,2]/A[2,2]))

###############################################################################
# Switching to multispecies case
# create the differential equation system
n_spec <- 10 # number of species in the system

# Lotka-Volterra model
A <- matrix(runif(n_spec^2,0,1), n_spec, n_spec) # interaction matrix
diag(A) <- runif(n_spec, n_spec, n_spec+1) # to ensure coexistence in the example
mu <- runif(n_spec,1,2) # intrinsic growth rate
test_f <- function(N){
  return(mu - A%*%N)
}

# compute relevant parameters with software
pars <- NFD_model(test_f, n_spec, from_R = TRUE)
ND <- pars$ND
NO <- pars$NO
FD <- pars$FD
c <- pars$c

# manualy check results for the two species case
# see appendix for proof of correctness
NO_check_m <- rep(0, n_spec)
FD_check_m <- rep(0, n_spec)
for (i in 1:n_spec){
  denominator = 0
  numerator = 0
  for (j in 1:n_spec){
    if (i==j) next
    numerator = numerator + pars$N_star[i,j]*A[i,j]
    denominator = denominator + pars$N_star[i,j]*sqrt(A[i,j]/A[j,i]*A[i,i]*A[j,j])
    NO_check_m[i] = numerator/denominator
    FD_check_m[i] = 1-denominator/mu[i]
  }
}
ND_check_m = 1-NO_check_m
###############################################################################
# passing additional arguments

n_spec <- 3 # number of species in the system

# Lotka-Volterra model
A <- matrix(runif(n_spec^2,0,1), n_spec, n_spec) # interaction matrix
diag(A) <- runif(n_spec, n_spec, n_spec+1) # to ensure coexistence in the example
mu <- runif(n_spec,1,2) # intrinsic growth rate
# LV model now depends on additional parameters
test_f <- function(N, mu, A){
  return(mu - A%*%N)
}

# we can pass additional arguments to the test function via the args
# argument. Pass an unnamed list (i.e. not list(mu=mu, A = A)), as
# python will not handle them correctly
pars <- NFD_model(test_f, n_spec, args = list(mu, A), from_R = TRUE)

###############################################################################
# passing additional starting estimates
# NFD_model computes automatically ND and FD for any model, however it's based
# on numerical solving, which requires starting estimates.
# the default for starting estimates is 1, so if your equilibrium abundance
# is very high, it will not find it.
n_spec <- 3
mu <- runif(n_spec, 1,2)*1e10

# will result in an error, as we don't help NFD_model
pars <- NFD_model(test_f, n_spec, args = list(mu, A), from_R = TRUE)

# we don't have to pass the equilibrium densities correctly, just valid
# starting estimates, e.g. N_star = 1e10 will suffice in this case
pars <- list(N_star = matrix(1e10, n_spec, n_spec))
pars <- NFD_model(test_f, n_spec, args = list(mu, A),pars = pars,
                  from_R = TRUE)

# N_star and c can be passed as initial starting guesses for c and N_star
# N_star and c must be matrices with shape (n_spec, n_spec)
# We can pass more specific starting guesses for each species absent
# N_star[1,] should be the densities of the community when species 1 is absent
# The value of N_star[1,1] will be ignored, no matter what
