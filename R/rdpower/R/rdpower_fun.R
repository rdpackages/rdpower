###################################################################
# Auxiliary functions for rdpower
# !version 2.2 20-Jun-2022
# Authors: Matias Cattaneo, Rocio Titiunik, Gonzalo Vazquez-Bare
###################################################################


#################################################################
# Power function
#################################################################

rdpower.powerfun <- function(n,tau,stilde,z){

  x <- 1 - pnorm(sqrt(n)*tau/stilde+z) + pnorm(sqrt(n)*tau/stilde-z)

  return(x)

}

#################################################################
# Power function derivative w.r.t. n
#################################################################

rdpower.powerfun.dot <- function(n,tau,stilde,z){

  x <- (dnorm(sqrt(n)*tau/stilde-z)-dnorm(sqrt(n)*tau/stilde+z))*tau/(2*stilde*sqrt(n))

  return(x)

}

#################################################################
# Power function derivative w.r.t. tau
#################################################################

rdpower.powerfun.dot.tau <- function(n,tau,stilde,z){

  x <- (dnorm(sqrt(n)*tau/stilde-z)-dnorm(sqrt(n)*tau/stilde+z))*sqrt(n)/stilde

  return(x)

}

#################################################################
# Newton-Raphson for sample size calculation
#################################################################

rdpower.powerNR <- function(x0,tau,stilde,z,beta){

  tol <- 1
  iter <- 0

  while (tol>.Machine$double.eps){
    iter <- iter + 1
    k <- 1
    # Check if derivative at x0 is too small
    while (rdpower.powerfun.dot(x0,tau,stilde,z)<.00001){
      x0 <- 1.2*x0*(rdpower.powerfun(x0,tau,stilde,z)<=beta) + .8*x0*(rdpower.powerfun(x0,tau,stilde,z)>beta)
      iter <- iter + 1
    }

    x1 <- x0 - (rdpower.powerfun(x0,tau,stilde,z)-beta)/rdpower.powerfun.dot(x0,tau,stilde,z)

    # Check if x1 is negative or too small
    while (x1<2){
      x1 <- x0 - k*(rdpower.powerfun(x0,tau,stilde,z)-beta)/rdpower.powerfun.dot(x0,tau,stilde,z)
      k <- k/2
      iter <- iter + 1
    }

    tol <- abs(rdpower.powerfun(x1,tau,stilde,z)-beta)
    x0 <- x1

  }

  b <- rdpower.powerfun(x1,tau,stilde,z)
  m <- ceiling(x1)

  output <- list(m = m,
                iter = iter,
                powercheck = b)

  return(output)
}

#################################################################
# Newton-Raphson for MDE calculation
#################################################################

rdpower.powerNR.mde <- function(n,tau0,stilde,z,beta){

  tol <- 1
  iter <- 0

  while (tol>.Machine$double.eps){
    iter <- iter + 1

    # Check if derivative at x0 is too small
    while (rdpower.powerfun.dot.tau(n,tau0,stilde,z)<.00001){
      x0 <- 1.5*tau0*(rdpower.powerfun(n,tau0,stilde,z)<=beta) + .5*tau0*(rdpower.powerfun(n,tau0,stilde,z)>beta)
      iter <- iter + 1
    }

    tau1 <- tau0 - (rdpower.powerfun(n,tau0,stilde,z)-beta)/rdpower.powerfun.dot.tau(n,tau0,stilde,z)


    tol <- abs(rdpower.powerfun(n,tau1,stilde,z)-beta)
    tau0 <- tau1

  }

  b <- rdpower.powerfun(n,tau1,stilde,z)

  output <- list(mde = tau1,
                iter = iter,
                powercheck = b)

  return(output)
}
