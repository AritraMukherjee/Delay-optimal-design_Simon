#############################################################################################

# Functions:delay_optimal()
# Last modified: 15 Sept 2021

############################################################################################

# Inputs

#    p0 - Success probability under the null hypothesis
#    p1 - Minimum success probability under the alternate hypothesis
# alpha - Desired type I error
#  beta - Desired type II error
#     t - Total recruitment time
#    m0 - Time to observe the primary outcome

#############################################################################################

# Output

# A vector of the values for (n1, n2, r1, r) for a delay optimal design assuming uniform and
# linear recruitments respectively along with their ESS

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

overrun_linear  <- function(n, n1, t, m0)
  {
  # Find delta
  delta <- 2*n/(t*(t+1))
  # Find t1
  t1    <- -1/2 + 1/2*sqrt(1 + 4*n1*t*(t + 1)/n)
  # Find the number of patients recruited
  n0    <- m0*delta*t1 + delta*m0*(m0 + 1)/2
  if (n0 <= n - n1) {
    n0
  } else {
    n - n1
  }
}

delay_optimal   <- function(p0, p1, alpha, beta, t, m0) {
  # From Fleming (1982)
  exact_n                <-
    ((stats::qnorm(1 - beta)*sqrt(p1*(1- p1)) +
        stats::qnorm(1 - alpha)*sqrt(p0*(1- p0)))/(p1 - p0))^2
  feasible_designs       <- saGS(2, p0, p1, alpha, beta, 1,
                                 ceiling(1.5*exact_n), 1, 0, 0, 0, 1)
  # ESS considering the delay
  # ESS_delay = (n1 + overrun)*(Prob that the trial stopped at first stage) +
  #             n*(Prob that the trial stopped at second stage)
  #           = ESS + overrun*PET
  ESS_uniform            <- ESS_linear <- numeric(nrow(feasible_designs))
  for(i in 1:nrow(feasible_designs)) {
    ESS_uniform[i]       <- feasible_designs[i, 10] +
      overrun_uniform((feasible_designs[i, 2] + feasible_designs[i, 1]),
                      feasible_designs[i, 1], t, m0)*feasible_designs[i, 8]
    ESS_linear[i]        <- feasible_designs[i, 10] +
      overrun_linear((feasible_designs[i, 2] + feasible_designs[i, 1]),
                     feasible_designs[i, 1], t, m0)*feasible_designs[i, 8]
  }
  feasible_designs       <- cbind(feasible_designs, ESS_uniform, ESS_linear)
  min_uniform            <- min(feasible_designs[, 16])
  argmin_uniform         <- logical(nrow(feasible_designs))
  for (i in 1:nrow(feasible_designs))
    {
    argmin_uniform[i]    <- isTRUE(all.equal(ESS_uniform[i], min_uniform))
    }
  if (sum(argmin_uniform) == 1)
    {
    null_optimal_uniform <- feasible_designs[argmin_uniform, ]
    }
  else
    {
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
  data                   <- c(null_optimal_uniform[1:4],
                              null_optimal_linear[1:4],
                              null_optimal_uniform[16], null_optimal_linear[17])
  names(data)            <- c("n1 (uniform)", "n2 (uniform)", "r1 (uniform)",
                              "r (uniform)", "n1 (linear)", "n2 (linear)",
                              "r1 (linear)", "r (linear)", "ESS (uniform)",
                              "ESS (linear)")
  data
}

