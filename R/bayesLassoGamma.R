#' @title bayesLassoGamma
#'
#' @description
#'
#' @param x Matrix of predictors, size n by p
#'
#' @param y Matrix of response, size p by 1
#'
#' @param Total Total number of iterations for Gibbs sampler
#'
#' @param B Burn-in iterations for Gibbs sampler
#'
#' @return Returns a matrix of beta coefficients from Gibbs sampler runs, size Total by p
#' @examples
#' x <- matrix(0, nrow=60, ncol=7)
#' e <- numeric(60)
#' sigma <- 2.5
#' for (i in 1:60){
#'  x[i,] <- stats::rnorm(7, 0, 1)
#'  e[i] <- stats::rnorm(1,0,sd=sigma^2)
#' }
#' y <- numeric(60)
#' for (i in 1:60){
#'  y[i] <- 3*x[i,4]+4*x[i,5]+e[i]
#' }
#' y <- y - mean(y) # centered y
#' bayesLassoGamma(x,y) # we expect only fourth and fifth coefficient to be far from zero

#' @export bayesLassoGamma
#' @importFrom
#' @references

bayesLassoGamma <- function(x,y,Total=10000, B=5000){
  n <- nrow(x)
  p <- ncol(x)

  # initial values
  xtx <- t(x) %*% x
  xy <- t(x) %*% y
  beta <- solve(xtx, xy) # OLS estimates
  tau_sq_inv <- 1 / beta^2 # initialize tau_sq_inv
  res <- drop(y - x %*% beta)
  sigma_sq <- mean(res^2) # initialize sigma_sq = MSE
  lambda <- p * sqrt(sigma_sq) / sum(abs(beta))
  r <-1
  theta <- 1.78

  # to store results
  beta_store <- matrix(0, nrow=Total, ncol=p)
  lambda_store <- numeric(Total)

  # Gibbs sampler
  for (i in 1:Total){
    # sample beta
    Dinv <- diag(as.vector(tau_sq_inv))
    Ainv <- solve(xtx + Dinv)
    mean <- Ainv %*% t(x) %*% y
    var <- sigma_sq*Ainv
    beta_star <- mnormt::rmnorm(n=1,mean, var) #library(mnormt)

    # sample sigma^2
    shape <- (n-1)/2+p/2
    scale <- (crossprod(y-x%*%beta_star))/2 + (t(beta_star) %*% Dinv %*% beta_star)/2
    sigma_sq_star <- 1/stats::rgamma(1, shape, scale)

    # sample 1/tau^2
    param1 <- sqrt(lambda^2*sigma_sq_star/beta_star^2)
    param2 <- lambda^2
    tau_sq_inv_star <- numeric(p)
    tau_sq_inv_star <- inversegaussian(param1, param2)

    # sample lambda: gamma prior
    shape <- p + r
    par <- sum(1/tau_sq_inv_star)/2 + theta
    lambda_star <- sqrt( stats::rgamma(1, shape, 1/par) )

    # update values
    beta <- beta_star
    sigma_sq <- sigma_sq_star
    tau_sq_inv <- tau_sq_inv_star
    lambda <- lambda_star

    # store results
    beta_store[i,] <- as.vector(beta)
    lambda_store[i] <- lambda
  }
  return(beta_store)
}

