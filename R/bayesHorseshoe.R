#' @title Function to perform Bayesian horseshoe regularization in linear regression
#'
#' @description Implementation of bhs MATLAB function.
#'
#' @param x Matrix of predictors, size n by p. Should be centered and scaled before calling the function.
#'
#' @param y Matrix of response, size p by 1. Should be centered before calling the function.
#'
#' @param Total Total number of iterations for Gibbs sampler that should be stored after thinning
#'
#' @param B Burn-in iterations for Gibbs sampler. Here burn-in iterations are not counted as part of Total iterations.
#'
#' @param thin Thinning parameter
#'
#' @return Returns a matrix of beta coefficients after discarding for burn-in iterations, size Total*p.
#' @examples
#' @export bayesHorseshoe
#' @references E. Makalic and D. F. Schmidt. (2016) A Simple Sampler for the Horseshoe Estimator.
#' IEEE Signal Processing Letters, 23(1), 179-182. Jan. 2016, doi: 10.1109/LSP.2015.2503725
#'
#' Rue, H. (2001). Fast Sampling of Gaussian Markov Random Fields. Journal of the Royal Statistical Society.
#' Series B (Statistical Methodology), 63(2), 325–338. http://www.jstor.org/stable/2680602
#'
#' Bhattacharya A., Antik Chakraborty A. & Mallick B. K. (2016). Fast sampling with Gaussian scale mixture priors in high-dimensional regression.
#' Biometrika, 103(4), 985–991. https://doi.org/10.1093/bi omet/asw042

bayesHorseshoe <- function (x,y, T=10000, B=5000, thin=11) {
  n <- nrow(x)
  p <- ncol(x)

  # initialize
  sigma_sq <- 1
  lambda_sq <- runif(p)
  lambda <- sqrt(lambda_sq)
  tau <- 1
  tau_sq <- 1
  epsilon <- 1 #xi
  v <- rep(1,p) #nu

  # to store
  beta_store <- matrix(NA, nrow = T, ncol = p )

  k <- 0
  iter <- 0
  # Gibbs sampler
  while(k<T){
    # sample beta (Rue 2001)   this also works !!!! this is same as in horsehoe apackage
    sigma <- sqrt(sigma_sq)
    lambda_star <- tau * lambda
    L <- chol((1/sigma_sq) * (xtx + diag(1/as.numeric(lambda_star^2), p, p)))
    a <- solve(t(L), t(t(y) %*% x)/sigma_sq)
    mu <- solve(L, a)
    u <- solve(L, rnorm(p))
    beta <- mu + u

    # sample sigma_sq
    param1 <- (n+p)/2
    e <- y- x %*% beta
    param2 <- ( t(e) %*% (e) )/2 + sum(((beta^2)/lambda_sq))/(tau_sq*2)
    sigma_sq <- 1/rgamma(1,param1,rate=param2)

    # sample lambda_sq
    value<- (1/v)+(beta)^2/(2*sigma_sq*tau_sq)
    lambda_sq <- 1/ rexp(p,rate=value)
    lambda <- sqrt(lambda_sq)

    # sample tau_sq
    param1 <- (p+1)/2
    param2 <- (1/epsilon)+(1/2)*sum(beta^2/lambda_sq)/sigma_sq
    tau_sq <- 1/rgamma(1,param1, rate=param2)
    tau <- sqrt(tau_sq)

    # sample v
    v <- 1/rexp(p,rate=1+(1/lambda_sq))

    # sample epsilon
    epsilon <- 1/rexp(1,rate=1+(1/tau_sq))

    #store results
    iter <- iter+1
    if (iter>B ){
      if(iter%%thin==0){
        k <- k+1
        beta_store[k,] <- as.vector(beta)
      }
    }
  }
  return(beta_store)
}
