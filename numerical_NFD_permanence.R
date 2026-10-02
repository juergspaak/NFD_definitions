## Numerical NFD in R
#
# @author: J.W.Spaak
# Numerically compute ND and FD for a model
# 

library(nleqslv)
library(numDeriv)
library(deSolve)

#### Main function: NFD_model ####

NFD_model_upgrade = function(f, n_spec = 2, args = c(), monotone_f = TRUE,
                             pars = NULL, f_jacobian = NULL,
                             xtol = 1e-8, equi_time = seq(1000), dN_dt = 1e-2,
                             N_min = 1e-3, full_equi = NULL,
                             plot = F, plot_dynamics = F,
                             pres = NULL){
  ##### Description ####
  #   Compute the ND and FD for a differential equation f
  #   
  #   Compute the niche difference (ND), niche overlap (NO), 
  #   fitness difference(FD) and conversion factors (c)
  #   
  #   Parameters
  #   -----------
  #   f : callable ``f(N, *args)``
  #       Percapita growth rate of the species.
  #       1/N dN/dt = f(N)
  #       
  #   n_spec : int, optional, default = 2
  #       number of species in the system
  #   args : tuple, optional
  #       Any extra arguments to `f`
  #   monotone_f : boolean or array of booleans (length: n_spec), default = True
  #       Whether ``f_i(N_i,0)`` is monotonely decreasing in ``N_i``
  #       Can be specified for each function separately by passing an array.
  #   pars : dict, default {}
  #       A dictionary to pass arguments to help numerical solvers.
  #       The entries of this dictionary might be changed during the computation
  #       
  #       ``N_star`` : ndarray (shape = (n_spec, n_spec))
  #           N_star[i] starting guess for equilibrium density with species `i`
  #           absent. N_star[i,i] is set to 0 
  #       ``r_i`` : ndarray (shape = n_spec)
  #           invsaion growth rates of the species
  #       ``c`` : ndarray (shape = (n_spec, n_spec))
  #           Starting guess for the conversion factors from one species to the
  #           other. `c` is assumed to be symmetric an only the uper triangular
  #           values are relevant
  #   xtol: float, default 1e-10
  #       Precision requirement of solving
  #   equi_time: vector, default seq(1000)
  #       Time for which ode maximally tries to find equilibrium
  #       Increase this if equilibrium is not found automatically
  #    dN_dt: float, default 1e-2
  #       Precision after which ode solver thinks equilibrium is found
  #       Decrease if equilibrium is not found automatically
  #    N_min, float, default 1e-3
  #        Exctintion threshhold, species below this density are assumed to be
  #        extinct. Change if you expect very low or very high equilibrium densities
  #    full_equi, array of float or NULL, default NULL
  #        Initial guess for equilibrium of full community.
  #    plot_dynamics, bool, default = FALSE
  #        If TRUE, then dynamics of all equilibria will be ploted.
  #    plot, boolean, default = FALSE
  #        If set to TRUE, will plot the resulting niche and fitness differences
  #   Returns
  #   -------
  #   pars : dict
  #       A dictionary with the following keys: 
  #           
  #   ``N_star`` : ndarray (shape = (n_spec, n_spec))
  #       N_star[i] equilibrium density with species `i`
  #       absent. N_star[i,i] is 0
  #   ``r_i`` : ndarray (shape = n_spec)
  #       invasion growth rates of the species
  #   ``c`` : ndarray (shape = (n_spec, n_spec))
  #       The conversion factors from one species to the
  #       other. 
  #   ``ND`` : ndarray (shape = n_spec)
  #       Niche difference of the species to the other species
  #       ND = (r_i - eta)/(mu -eta)
  #   ``FD`` : ndarray (shape = n_spec)
  #       Fitness difference according to Spaak et al. 2021
  #       FD = (0 - eta)/(mu - eta)
  #   ``eta``: ndarray (shape = n_spec)
  #       no-niche growth rate f(\sum c_j^i N_j^(-i),0)
  #       eta and fc are identical, but both are maintained for compatibility
  #   ``mu``: ndarray (shape = n_spec)
  #       intrinsic growth rate f(0,..., 0)
  #   ``surv``: ndarray (shape = n_spec, type = bool)
  #       Boolean array which is TRUE for surviving species and F for absent species
  #   
  #   Raises:
  #       InputError:
  #           Is raised if system cannot automatically solve equations.
  #           Starting estimates for N_star and c should be passed.
  #   
  #   Examples:
  #       See "Example,compute NFD.R"
  #   
  #   Debugging:
  #       If an Error is raised the problem causing information is saved in
  #       pars.
  #       To access it rerun the code in the following way (or similar)
  #           
  #       pars = new.env()
  #       pars = NFD_model(f, pars = pars)
  #       print(as.list(pars))
  #       
  #       pars will then contain additional information
  #       
  #   Literature:
  #   "Intuitive and broadly applicable definitions of 
  # niche and fitness differences", J.W.Spaak, F. deLaender
  #   DOI: https://doi.org/10.1101/482703 
  #   
  
  ##### Code ####
  if (n_spec == 1){
    # species case, single species are assumed to have ND = 1
    stop(
      paste("ND and FD are not (properly) defined for a single",
            "species community.",
            "If needed assign manualy ND = 1 and FD = 0 for this case", 
            sep = "\n"))
  }
  
  equi_args <- list(time = equi_time,
                    dN_dt = dN_dt,
                    N_min = N_min,
                    xtol = xtol,
                    plot_dynamics = plot_dynamics,
                    full_equi = full_equi)
  # check whether all input is parsed correctly
  equi_args <- input_check(n_spec, f, args, equi_args, pars)
  
  # prepare pars list and check its input
  pars <- preconditioner_upgrade(f, args, n_spec, pars, xtol, monotone_f, pres)
  
  # find equilibrium of the entire community
  equi <- find_equi(equi_args$full_equi[pars$pres], pars$pres, pars, equi_args)
  
  # did we find the equilibrium
  surv <- equi$equi != 0 # all surviving species
  
  # did at least one species survive?
  if (!any(surv)) stop("Apparently no species survived, please check your model",
                       ". Potentially rerun the command with `plot_dynamics = T`.")
  
  # store which species survive
  pars$surv <- surv # which species survive, initially assume all
  pars$equi <- equi$equi
  
  # invasion community for all non-surviving species is the full community
  if (any(!pars$surv)){
    
    pars$N_star[!pars$surv,] = matrix(pars$equi,
                                      sum(!pars$surv), n_spec, byrow =T)
  }
  ind_surv <- seq(n_spec)[pars$surv] # indices of surviving species

  # find other equilibria
  for (i in ind_surv){
    # find the resident community of the absent species
    pres_temp <- pars$surv
    pres_temp[i] = F # remove focal species
    if (any(pres_temp)){
      equi_temp <- find_equi(pars$N_star[i,pres_temp], pres_temp, pars, equi_args)
      pars$N_star[i,] = equi_temp$equi
    }
    else pars$N_star[i,] <- 0
  }
  
  # for each species compute the invasion growth rate
  for (i in 1:n_spec){
    pars$r_i[i] = pars$f(pars$N_star[i,])[i]
  }
  if (any((pars$r_i[pres]>0) != pars$surv[pres])) stop(paste0(
    "Automatically computed equilibria lead to wrong coexistence.",
    "To solve this problem, please provide equilibrium estimates ",
    "via `pars$N_star` argument. Alternatively, please increase equi_time, ",
    "decrease N_min, and/or decrease dN_dt."))
  
  if (any((pars$N_star >0) & (pars$N_star<N_min))){
    warning(paste0("Some species have equilibrium density below proposed minimal density ",
                   "`N_min`. This may be fine but please double check.",
                   "Decreasing `N_min` will get rid of this warning."))
  }
  # list of all species
  l_spec <- 1:n_spec
  
  for (i in 1:n_spec) {
    for (j in 1:n_spec) {
      if (i >= j) {  # c is assumed to be symmetric, c[i,i] = 1
        next  # Equivalent to 'continue' in Python
      }
      conversion <- solve_c(pars, c(i, j), monotone_f, xtol = xtol)
      pars$c[i,j] <- conversion[1]
      pars$c[j,i] <- conversion[2]
    }
  }
  
  ## Compute eta, ND and FD
  pars$eta <- numeric(n_spec)
  for (i in l_spec) {
    sp <- c(i, l_spec[l_spec != i])
    pars$eta[i] <- pars$f(switch_niche(pars$N_star[i,], sp,
                                       pars$c[i,-i]))[i]
  }
  
  # prepare returning values
  pars$ND = (pars$r_i - pars$eta)/(pars$mu - pars$eta)
  pars$FD = (0-pars$eta)/(pars$mu - pars$eta)
  
  # assing niche and fitness differences to species that don't survive
  pars$ND[pars$eta == pars$mu] = 1
  pars$FD[pars$eta == pars$mu] = 1
  
  
  if (plot) plot_NFD(pars)
  
  return(pars)
  
  
}

#### Supplementary functions ####

####### identical to other code:
input_check = function(n_spec, f, args, equi_args, pars) {
  # check input on (semantical) correctness
  if (!(round(n_spec) == n_spec)) {
    stop("Number of species (`n_spec`) must be an integer")
  }

  # check whether `f` is a function
  tryCatch({
    mu <- do.call(f, c(list(rep(0, n_spec)), args))
    if (length(mu) != n_spec) {
      stop("`f` must return an array of length `n_spec`")
    }
  }, error = function(e) {
    message("function call of `f` did not work properly")
    stop(e)
  })
  
  if (any(!is.finite(mu))) {
    stop("For some species f(0) seems to be not defined.")
  }
  
  # check equi_args must be positive real numbers
  for (key in c("xtol", "dN_dt", "N_min")){ # 
    p <- equi_args[[key]]
    if (!is.numeric(p)
        || length(p) != 1
        || !is.finite(p)
        || Im(p) != 0
        || Re(p) <= 0) {
      stop(sprintf("%s must be a single positive real number", key))
    }
  }
  
  # plot_dynamics must be a boolean
  if (!is.logical(equi_args$plot_dynamics)
      || length(equi_args$plot_dynamics) !=1){
    stop("`plot_dynamics must be TRUE or FALSE")
  }
  
  # equi_time should be a numerical vector or timepoint
  if (!is.numeric(equi_args$time)){
    stop("`equi_time` must be numerical")
  }
  if (length(equi_args$time) == 1){
    equi_args$time <- seq(0, abs(equi_args$time),
                          length.out = 101) # convert to sequence
  }
  
  # full equi must be positive and real numbers
  if (is.null(equi_args$full_equi)) equi_args$full_equi <- rep(1, n_spec)
  else if (is.numeric(equi_args$full_equi) # is numeric matrix
           && length(equi_args$full_equi) # correct length
           && all(Im(equi_args$full_equi) == 0) # is a positive real number
           && all(equi_args$full_equi>=0)
           && all(is.finite(equi_args$full_equi))) {} # everything fine
  else stop(paste0("`full_equi` must be the starting density of the",
                   " simulation, i.e. a vector of length `n_spec` of", 
                   " positive real numbers."))
  
  
  # returns nothing
  return(equi_args)
}

##### preconditioner() function #####
preconditioner_upgrade = function(f, args, n_spec, pars, xtol, monotone_f, pres){
  ###### Description ######
  
  # Returns equilibria densities and invasion growth rates for system `f`
  #   
  #   Parameters
  #   -----------
  #   same as `find_NFD`
  #           
  #   Returns
  #   -------
  #   pars : dict
  #       A dictionary with the keys:
  #       
  #       ``N_star`` : ndarray (shape = (n_spec, n_spec))
  #           N_star[i] is the equilibrium density of the system with species 
  #           i absent. The density of species i is set to 0.
  #       ``r_i`` : ndarray (shape = n_spec)
  #           invsaion growth rates of the species
  #   
  
  
  if(is.null(pars)){
    pars = list()
  }
  pars$n_spec <- n_spec # species richness
  pars$xtol <- xtol # tolerance
  
  # assume all species are present if not specified otherwise
  if (is.null(pres)) pars$pres <- rep(T, n_spec)
  else pars$pres <- pres
  
  # check whether pres is meaningful
  if (!is.logical(pars$pres)
      || length(pars$pres) !=n_spec){
    stop("`pres` must be a vector of length `n_spec` with TRUE and FALSE")
  }
  if (!any(pars$pres)) stop("At least one species must be present in `pres`")
  
  
  warn_string <- "pars[%s] must be array with shape %s and real positive numbers. The values will be computed automatically"
  
  # check given keys of pars for correctness
  for (key in c("c", "N_star")) {
    
    if (is.null(pars[[key]])) {
      # key not present in `pars`
      pars[[key]] <- matrix(1, n_spec, n_spec)
    } else if (is.matrix(pars[[key]]) # input is a matrix
               && all(dim(pars[[key]]) == c(n_spec, n_spec)) # of the right shape
               && is.numeric(pars[[key]]) # is numeric matrix
               && all(Im(pars[[key]]) == 0) # is a positive real number
               && all(pars[[key]]>=0)
               && all(is.finite(pars[[key]]))) pars[[key]] # values are fine
    else{
      pars[[key]] <- matrix(1, n_spec, n_spec) # replace with ones
      # warn user
      warning(sprintf(warn_string, key, paste(c(n_spec, n_spec), collapse = " ")))
    }  
  }
  
  pars$f <- function(N) {
    # Allow passing infinite species densities to per capita growth rate
    if (any(is.infinite(N))) {
      return(rep(-Inf, length(N)))
    } else {
      N <- N
      N[N < 0] <- 0 # Function might be undefined for negative densities
      growth  = do.call(f, c(list(N = N), args))
      if(!(all(is.finite(growth)))){
        pars$N_error <- N
        pars$f_error <- growth
        stop(paste0("Function returns non-finite values at densities: (",
                    paste(N, collapse = ", "),
                    ")"))
      }
      return(do.call(f, c(list(N = N), args)))
    }
  }
  
  # Monoculture growth rate
  pars$mu <- pars$f(rep(0, n_spec))
  
  return(pars)
}

find_equi <- function(N_start, present, pars, equi_args){
  # which community are we working with?
  com = ifelse(all(present == pars$pres), "full-",
               paste0("-", (1:pars$n_spec)[present != pars$surv], " "))

  # when to stop numerical integration
  stop_ode <- function(t, N, parms){
    N_all <- rep(0, pars$n_spec) # add densities of absent species
    N_all[present] = N
    stop_cond <- (abs(pars$f(N_all))<equi_args$dN_dt |# is species close to equilibrium?
                    (N_all<equi_args$N_min & # or low density and negative growth rate
                       pars$f(N_all)<0) |
                    (!present)) # or simply absent
    # returns 0 if all species are at equilibrium        
    1 - min(stop_cond)
  }
  
  # convert per-capita growth rate to actual growth rate
  ode_fun <- function(t, N, nothing ){
    N_all <- rep(0, pars$n_spec) # add densities of absent species
    N_all[present] = N
    # use actual growth rates, not per-capita
    return(list(pars$f(N_all)[present]*N)) # only return present species
  }
  
  if (!stop_ode(0, N_start, nothing)){
    # initial estimate already close to equilibrium
    N_pre = rep(0, pars$n_spec)
    N_pre[present] = N_start
  }
  else{
    out <- ode(N_start, equi_args$time, ode_fun, parms = list(),
               rootfun = stop_ode)
    # fill in approximate equilibrium
    N_pre = rep(0, pars$n_spec)
    N_pre[present] = out[nrow(out), 2:ncol(out)]
  }
  
  # remove all absent species from equilibrium
  present_new = N_pre>equi_args$N_min
  N_pre[!present_new] <- 0 # remove species that went extinct
  fsolve_fun <- function(N){
    N_all <- rep(0, pars$n_spec) # add densities of absent species
    N_all[present_new] = N
    # use actual growth rates, not per-capita
    return(N*pars$f(N_all)[present_new]) # only return present species
  }
  
  # did at least one species survive?
  if (!any(present_new)){
    # check whether all species have negative intrinsic growth rate
    
    warning(paste0("All species went extinct for ", com, " community."))
    if (equi_args$plot_dynamics && exists("out")){
        plot_equi(out, N_pre[present], present_new[present], equi_args, com)}
    
    if (any(pars$f(N_pre)[present]>equi_args$xtol)){
      if (exists("out")){
        plot_equi(out, N_pre[present], present_new[present], equi_args, com)}
      else {
        warning(paste0("Can't plot dynamics because species already",
                       " at equilibrium for ", com, " community."))
      }
      stop(
      paste0("Could not find automatical equilibrium for ", com, "community. ",
           "Please check plot. To solve this problem, please provide ",
           "equilibrium estimates via `pars$N_star` or `full_equi` argument. ",
           "Alternatively, please increase equi_time, decrease N_min,
                and/or decrease dN_dt."))}
    
    return(list(equi = N_pre, equi_found = T))
  }
  
  # check result using root finder close to end of ODE function
  tryCatch({fsolve_result <- nleqslv(N_pre[present_new], fsolve_fun,
                           control = list(ftol = equi_args$xtol,
                                          xtol = equi_args$xtol))
           N_pre[present_new] <- fsolve_result$x},
    error = function(e){}) # can't update equilibrium
  
  # Calculate the Jacobian manually using numerical differentiation
  # jacobian of the differential equation (not the per-capita growth rates)
  jac <- numDeriv::jacobian(fsolve_fun, N_pre[present_new])
  
  # species that were present in simulation, but went extinct
  extinct = !present_new # all species absent in N_pre
  # species that were initially not present are removed from extinct
  extinct[!present] = F 
  equi_found <- (max(abs(pars$f(N_pre)[present_new])<equi_args$xtol) # is an equilibrium
                 & all(N_pre[present_new] > 0) # is feasible
                 & all(is.finite(N_pre)) # isn't nonsense
                 & max(Re(eigen(jac)$values)) < 0 # is stable
                 & all(pars$f(N_pre)[extinct] < 0) # absent species can't invade
  )
  
  if (equi_args$plot_dynamics
      || (!equi_found)){
    if (exists("out")){
      plot_equi(out, N_pre[present], present_new[present], equi_args, com)}
    else {
      warning(paste0("Can't plot dynamics because species already",
                     " at equilibrium for ", com, " community."))
    }}
  
  if (!(equi_found)){
    stop(paste0("Could not find automatical equilibrium for ", com, "community. ",
                "Please check plot. To solve this problem, please provide ",
                "equilibrium estimates via `pars$N_star` or `full_equi` argument. ",
                "Alternatively, please increase equi_time, decrease N_min,
                and/or decrease dN_dt."))
  }
  
  return(list(equi = N_pre, equi_found = equi_found, jac = jac,
              no_invasion = pars$f(N_pre)))
  
}

plot_equi <- function(out, equi, loc_surv, equi_args, com){
  plot(NA, xlim = range(out[,1]),
       ylim = c(equi_args$N_min/100, max(out[,2:ncol(out)])*1.5),
       log = "y", main = paste0(com, "community"), xlab = "Time", ylab = "Density")
  for (i in 2:ncol(out)){
    lines(out[,1], out[,i], col = ifelse(loc_surv[i-1], "black", "red"))
    points(max(out[,1]), equi[i-1],
           col = ifelse(loc_surv[i-1], "black", "red"))
  }
  abline(h = equi_args$N_min, col = "red")
}

##### solve_c() function #####
solve_c <- function(pars, sp = c(1, 2), monotone_f = TRUE, xtol = 1e-10) {
  
  ###### Description ######
  
  # Find the conversion factor c for species sp
  #   
  #   Parameters
  #   ----------
  #   pars : dict
  #       Containing the N_star and r_i values, see `preconditioner`
  #   sp: array-like
  #       The two species to convert into each other
  #       
  #   Returns
  #   -------
  #   c : float, the conversion factor c_sp[0]^sp[1]
  # 
  
  ###### Code ######
  
  # get the niche overlap functions
  NO_case1 <- NO_fun(pars, sp)
  NO_case2 <- NO_fun(pars, rev(sp))
  NO_funs <- list(NO_case1[[1]],
                  NO_case2[[1]])
  special_case <- c(NO_case1[[2]], NO_case2[[2]])
  
  pars$NO_funs <- NO_funs
  # Do species interact?
  if (any(special_case == "No interaction")) {
    return(special_case_nointeraction(special_case == "No interaction", sp))
  }
  
  # Has one species reached minimal growth rate?
  if (any(special_case == "Mortality reached + facilitation")) {
    return(special_case_mort(special_case == "Mortality reached + facilitation", sp))
  }
  
  # Define the function to be solved
  inter_fun <- function(conversion) {
    NO_ij <- abs(NO_funs[[1]](conversion))
    NO_ji <- abs(NO_funs[[2]](1 / conversion))
    return(NO_ij - NO_ji)
  }
  
  # Use a generic numerical solver when `f` is not monotone
  if (!monotone_f) {
    conversion <- nleqslv(pars$c[sp[1], sp[2]], inter_fun, control = list(xtol = xtol))$x
    if (abs(inter_fun(conversion)) > xtol) {
      pars$`c found by nleqslv` <- conversion
      stop(paste0("Not able to find c_", sp[1], "^", sp[2], ".",
                  " Please pass a better guess for c_i^j via the `pars` argument"))
    }
    return(c(conversion, 1/conversion))
  }
  
  # If `f` is monotone, then the solution is unique, find it with a more robust method
  
  # Find interval for uniroot method
  a <- pars$c[sp[1], sp[2]]
  direction <- sign(inter_fun(a))
  if (direction == 0) {
    return(c(a, 1/a)) # initial c_value is correct
  }
  
  fac <- 2^direction
  if (!is.finite(direction)) {
    pars$`function inputs` <- lapply(list(sp, rev(sp)), function(es) {
      sapply(c(0, a, 1 / a), function(conversion) switch_niche(pars$N_star[es[1]], es, conversion))
    })
    pars$`function outputs` <- sapply(pars$`function inputs`, function(inp) pars$f(inp))
    stop("Function `f` seems to be returning nonfinite values")
  }
  
  b <- as.numeric(a * fac)
  
  while (sign(inter_fun(b)) == direction) {
    a <- b
    b <- b * fac
    if (!((2 * a == b) || (2 * b == a)) || sign(b - a) != direction) {
      stop(paste0("Not able to find c_", sp[1], "^", sp[2], ".",
                  " Please pass a better guess for c_i^j via the `pars` argument",
                  ". Please also check for non-positive entries in pars$c"))
    }
  }
  
  # Solve the equation using uniroot
  conversion <- tryCatch({
    uniroot(inter_fun, interval = c(a, b), tol = xtol/100)$root
  }, error = function(e) {
    stop("Fun
         ction does not seem to be monotone. Please run with `monotone_f = FALSE`")
  })
  
  # Test whether c actually is correct
  if (inter_fun(conversion) > xtol) {
    warning(sprintf(paste0("The computed values for c_(%d,%d) are above the",
                           " tolerance given by the user. Consider adjusting",
                           " xtol or pass with pars$c argument"), sp[1], sp[2]))
  }
  return(c(conversion, 1/conversion)) # Return c_i and c_j = 1/c_i
}

##### switch_niche() function #####
switch_niche <- function(N, sp, conversion = 0) {
  # switch the niche of sp[2:length(sp)] into niche of sp[1]
  N <- N # copy N to avoid modifying the original data
  N[sp[1]] <- N[sp[1]] + sum(conversion * N[sp[-1]], na.rm = TRUE)
  N[sp[-1]] <- 0
  return(N)
}

##### NO_fun() function #####
NO_fun <- function(pars, sp) {
  # returns the niche overlap function as a function of conversion
  # first check whether species are present
  if (pars$N_star[sp[1], sp[2]]>0){
    # growth when species sp[2] is removed
    f0 <- pars$f(switch_niche(pars$N_star[sp[1],], sp))[sp[1]]
    
    # growth when species sp[2] is put into niche of sp[1]
    fc <- pars$f(switch_niche(pars$N_star[sp[1],], sp, c_ij))[sp[1]]
    
    # check for special cases
    # has species reached mortality rate?
    if ((abs(f0-pars$r_i[sp[1]])<pars$xtol)
        & (abs(f0-fc)<pars$xtol)){
      warning(sprintf("Species %d reached mortality/minimal growth rate. This may result in special behaviour, please check c, ND and FD.", sp[1]))
      return(list(function(c) 1, "Mortality reached"))
    }
    
    if (abs(f0-pars$r_i[sp[1]])<pars$xtol){
      warning(sprintf("Species %d does not seem to affect species %d. This may result in nonfinite c, ND and FD values.", sp[2], sp[1]))
      return(list(function(c) 0, "No interaction"))
    }
    
    if ((abs(f0-fc)<pars$xtol)){
      warning(sprintf("Species %d reached mortality/minimal growth rate and is facilitate by species $d. This may result in special behaviour, please check c, ND and FD.", sp[1], sp[2]))
      return(list(function(c) Inf, "Mortality reached + facilitation"))
    }
    # return a function that computes niche overlap
    NO_pre = function(c_ij)
    {fc <- pars$f(switch_niche(pars$N_star[sp[1],], sp, c_ij))[sp[1]]
    return(abs((f0 - pars$r_i[sp[1]])/(f0 - fc)))}
  }
  else{ # N_star[i,j] = 0, hence use jacobian
    if (is.null(pars$f_jacobian)){# use numerical jacobian
      side = rep(1, pars$n_spec) # avoid evaluating negative densities
      jac = numDeriv::jacobian(pars$f, pars$N_star[sp[1],], side = side)}
    else {
      # use analytical jacobian
      jac <- tryCatch(do.call(f_jacobian, c(list(N = pars$N_star[sp[1],]), args)),
                      error = function(e) {
                        warning(sprintf(paste0("An error occured while ",
                                               "evaluating the jacobian at the equilibrium density",
                                               "for species %d. The error is stored in pars$error.",
                                               "The jacobian is now computed automatically."), sp[1]))
                        pars$error <- e
                        # compute numerically
                        numDeriv::jacobian(pars$f, pars$N_star[sp[1],], side = NULL)
                      })
    }
    NO_pre <- function(c_ij) return(1/c_ij*jac[sp[1], sp[2]]/jac[sp[1], sp[1]])
  }
  
  # save version of the function being 1 if 0/0
  NO_final <- function(c_ij){
    NO_value <- NO_pre(c_ij)
    ifelse(is.nan(NO_value), 1, NO_value)
  }
  return(list(NO_final, "standard"))
}



##### special_case() function #####
special_case_nointeraction <- function(no_comp, sp) {
  # Return c for special case where one species is not affected by competition

  if (all(no_comp)) {
    return(c(0, 0))  # Species do not interact at all, c set to zero
  } else if (all(no_comp == c(TRUE, FALSE))) {
    return(c(0, Inf))  # Only first species affected
  } else if (all(no_comp == c(FALSE, TRUE))) {
    return(c(Inf, 0))  # Only second species affected
  }
}

##### special_case_mort() function #####
special_case_mort <- function(mort, sp) {
  # Return c for special case where one species is not affected by itself

  if (all(mort)) {
    return(c(1,1))  # Both species have reached mortality rate
  } else if (all(mort == c(TRUE, FALSE))) {
    return(c(Inf, 0))  # Only first species affected
  } else if (all(mort == c(FALSE, TRUE))) {
    return(c(0, Inf))  # Only second species affected
  }
}

plot_NFD <- function(pars){
  xlim = range(c(pars$ND, -0.1, 1.1))
  
  plot(pars$ND, pars$FD, xlim = range(c(pars$ND, -0.1, 1.1)),
       ylim = range(c(pars$FD, -0.1, 1.1)),
       xlab = "Niche differences", ylab = "Fitness differences",
       col = ifelse(pars$r_i>0, "blue", "red"),
       pch = ifelse(pars$pres, 17, 16))
  if (!all(pars$pres)) {legend("topleft", legend = c("absent", "present"),
         pch = c(16, 17), bty = "n")}
  abline(h = c(0,1))
  abline(v = c(0,1))
  abline(0, 1, col = "red")
}