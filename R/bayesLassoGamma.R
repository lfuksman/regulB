#' @title Function to perform Bayesian lasso in linear regression
#'
#' @description
#' Consider the model \eqn{y = x\beta+\epsilon} where \eqn{y = (y_1,...,y_n)} are \eqn{n} responses
#' and \eqn{\epsilon \sim N(0,\sigma^2)}.
#' The function is based on Park and Casella (2008) and computes \eqn{\beta} estimates from posterior distribution.
#' Gibbs sampler generates \eqn{\beta} from multivariate normal distribution,
#' \eqn{\sigma^2} from inverse gamma, \eqn{\tau_i^{-2}} from inverse gaussian and
#' \eqn{\lambda^2} from gamma distribution. The algorithm for sampling from inverse gaussian
#' distribution is based on Michael, Schucany and Haas (1976).
#'
#' @param x Matrix of predictors, size n by p. Should be centered and scaled before calling the function.
#'
#' @param y Matrix of response, size p by 1. Should be centered before calling the function.
#'
#' @param T Total number of iterations for Gibbs sampler
#'
#' @param B Burn-in iterations for Gibbs sampler.Here burn-in ietrations are not counted as part of Total iterations.
#'
#' @return Returns a matrix of beta coefficients after discarding burn-in iterations,
#' size Total*p.
#'
#' @details Gibbs sampler is based on the following full conditional distributions:
#'
#' \eqn{\beta|. \sim MN(A^{-1}X'y,\sigma^2 * A^{-1}) } where
#' \eqn{ A=X'X+{D}_\tau^{-1} }
#'
#' \eqn{ \sigma^2|. \sim IG( (1/2) * (n-1+p), (1/2) * [(y-X\beta)'(y-X\beta) + \beta'{D}_\tau^{-1}\beta] ) }
#'
#' \eqn{ \tau_i^{-2}|. \sim Inverse-Gaussian (sqrt(\lambda^2 * \sigma^2 / \beta_i^2), \lambda^2) }
#'
#' \eqn{\lambda^2|. \sim Gamma(p+r, \delta+\sum_{i=1}^{p}[ (\tau_i^2)/2) ] ) }
#'
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
#' @references
#' Park T, & Casella G. (2008). The Bayesian Lasso. Journal of the American Statistical Association,
#' 103:482, 681-686, DOI: 10.1198/016214508000000337
#'
#' Michael, John R.; Schucany, William R.; Haas, Roy W. (1976),
#' "Generating Random Variates Using Transformations with Multiple Roots",
#' The American Statistician, 30 (2): 88–90, doi:10.1080/00031305.1976.10479147, JSTOR 2683801

bayesLassoGamma <- function(x,y,T=10000, B=5000){
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
  beta_store <- matrix(0, nrow=(T+B), ncol=p)
  lambda_store <- numeric(T)

  # Gibbs sampler
  for (i in 1:(T+B)){
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
  return(beta_store[-(1:B),])
}


