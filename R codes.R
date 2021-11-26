#############################################################################################

# Functions:delay_optimal()
# Last modified: 15 Sept 2021

############################################################################################

# Inputs

#    p0 - Success probability under the null hypothesis
#    p1 - Minimum success probability under the alternate hypothesis
# alpha - Desired type I error
#  beta - Desired type II error
#     k - Proportion of the total recruitment time when there is a linear recruitement 
#     t - Total recruitment time
#    m0 - Time to observe the primary outcome


#############################################################################################

# Output

# A vector of the values for (n1, n2, r1, r) for a delay optimal design assuming uniform, linear
# and mixed recruitments respectively along with their ESS

############################################################################################

# Required packages

#install.packages("Rcpp")

############################################################################################

# First need to source simon_cpp_code.cpp (ensure this is in the current working directory)
Rcpp::sourceCpp("./simon_cpp_code.cpp")

# Function for computing the number of patients recruited during data accrual
# assuming uniform and linear recruitments

overrun_uniform <- function(n, n1, t, m0)
  {
  # Find lambda, the rate of arrival of patients under a Uniform recruitment
  lambda <- n/t
  # Find the number of patients recruited
  n0     <- lambda*m0
  if (n0 <= n - n1) {
    n0
  } else {
    n - n1
  }
}

overrun_mixed  <- function(n, n1, k, t, m0) 
{
  if(k>0)
  {
  #Let us assume, patients are recruited in a linearly increasing pattern (as δt) 
  #upto k times of the total recruitment length (t). [0<k≤1]
  
  #Then, upto time point kt, the total recruitment= δ [1+2+ … + kt] 
  #=δkt(kt+1)/2
  #If we assume that, thereafter for (1-k)t timepoints patients are recruited uniformly 
  #at a rate of δkt, then the total recruitment for that period = δkt(1-k)t
  
  # For a Simon’s design, then, the total recruitment n should be equal to this total recruitment.
  # i.e. δkt(kt+1)/2+δkt(1-k)t=n
  
  # Find delta
  delta <- 2*n/(k*t*(1+(2-k)*t))
  # Find t1
  t1    <- -1/2 + 1/2*sqrt(1 + 8*n1/delta)
  
  # Find the number of patients recruited
  # case 1. t1 still falls under linear recruitment
  if(t1<=k*t)
  {
    # case i. t1+m0<=kt then recruitment is linear 
    if(t1+m0<=k*t) 
    {
      n0=m0*delta*t1 + delta*m0*(m0 + 1)/2
    }
    
    # case ii. t1+m0>kt; after t1 units to kt units of time, the recruitment is linear 
    # and from kt to m0 it's uniform
    else 
    {
      n0=delta*((k*t-t1)*t1+1/2*(k*t-t1)*(k*t-t1+1))+delta*k*t*(t1+m0-k*t)
    }
  }
  #case 2. at t1 uniform recruitment has started
  else
  {
    n0=m0*delta*k*t
  }
  
  if (n0 <= n - n1) 
    return(n0)
  else 
    return(n - n1)
  }
  
  if(k==0)
    n0=overrun_uniform(n,n1,t,m0)
}



delay_optimal   <- function(p0, p1, alpha, beta, k, t, m0) {
  # From Fleming (1982)
  exact_n                <-
    ((stats::qnorm(1 - beta)*sqrt(p1*(1- p1)) +
        stats::qnorm(1 - alpha)*sqrt(p0*(1- p0)))/(p1 - p0))^2
  feasible_designs       <- saGS(2, p0, p1, alpha, beta, 1,
                                 ceiling(1.5*exact_n), 1, 0, 0, 0, 1)
  #ESS considering the delay
  #ESS_delay = (n1 + overrun)*(Prob that the trial stopped at first stage) +
  # n*(Prob that the trial stopped at second stage)
  #   = ESS + overrun*PET
  ESS_uniform            <- ESS_linear     <- ESS_mixed   <- numeric(nrow(feasible_designs))
  for(i in 1:nrow(feasible_designs)) {
    ESS_uniform[i]       <- feasible_designs[i, 10] +
      overrun_uniform((feasible_designs[i, 2] + feasible_designs[i, 1]),
                      feasible_designs[i, 1], t, m0)*feasible_designs[i, 8]
    ESS_linear[i]        <- feasible_designs[i, 10] +
      overrun_mixed((feasible_designs[i, 2] + feasible_designs[i, 1]),
                     feasible_designs[i, 1],1, t, m0)*feasible_designs[i, 8]
    ESS_mixed[i]        <- feasible_designs[i, 10] +
      overrun_mixed((feasible_designs[i, 2] + feasible_designs[i, 1]),
                    feasible_designs[i, 1], k, t, m0)*feasible_designs[i, 8]
  }
  feasible_designs       <- cbind(feasible_designs, ESS_uniform, ESS_linear, ESS_mixed)
  min_uniform            <- min(feasible_designs[, 16])
  argmin_uniform         <- logical(nrow(feasible_designs))
  for (i in 1:nrow(feasible_designs)) {
    argmin_uniform[i]    <- isTRUE(all.equal(ESS_uniform[i], min_uniform))
  }
  if (sum(argmin_uniform) == 1) {
    null_optimal_uniform <- feasible_designs[argmin_uniform, ]
  } else {
    null_optimal_uniform <-
      feasible_designs[which(feasible_designs[argmin_uniform, 10] ==
                               min(feasible_designs[argmin_uniform, 10])), ]
  }
  min_linear             <- min(feasible_designs[, 17])
  argmin_linear          <- logical(nrow(feasible_designs))
  for (i in 1:nrow(feasible_designs)) {
    argmin_linear[i]     <- isTRUE(all.equal(ESS_linear[i], min_linear))
  }
  if (sum(argmin_linear) == 1) {
    null_optimal_linear  <- feasible_designs[argmin_linear, ]
  } else {
    null_optimal_linear  <-
      feasible_designs[which(feasible_designs[argmin_linear, 10] ==
                               min(feasible_designs[argmin_linear, 10])), ]
  }
  
  min_mixed             <- min(feasible_designs[, 18])
  argmin_mixed          <- logical(nrow(feasible_designs))
  for (i in 1:nrow(feasible_designs)) {
    argmin_mixed[i]     <- isTRUE(all.equal(ESS_mixed[i], min_mixed))
  }
  if (sum(argmin_mixed) == 1) {
    null_optimal_mixed  <- feasible_designs[argmin_mixed, ]
  } else {
    null_optimal_mixed  <-
      feasible_designs[which(feasible_designs[argmin_mixed, 10] ==
                               min(feasible_designs[argmin_mixed, 10])), ]
  }
  
  
  data                   <- c(null_optimal_uniform[1:4],
                              null_optimal_linear[1:4],
                              null_optimal_mixed[1:4],
                              null_optimal_uniform[16], null_optimal_linear[17],null_optimal_mixed[18])
  names(data)            <- c("n1 (optimal)", "n2 (uniform)", "r1 (uniform)",
                              "r (uniform)", "n1 (linear)", "n2 (linear)",
                              "r1 (linear)", "r (linear)", "n1 (mixed)", 
                              "n2 (mixed)", "r1 (mixed)", "r (mixed)","ESS (uniform)",
                              "ESS (linear)", "ESS(mixed)")
  data
}
