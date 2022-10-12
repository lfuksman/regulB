#' @title plotBayesGammaLasso
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
#' @return
#' @examples

#' @export plotBayesGammaLasso
#' @importFrom
#' @references

plotBayesGammaLasso <- function(x,y,Total=10000,B=5000){
  b <- bayesLassoGamma(x,y,Total,B)
  p <- ncol(x)
  b <- b[-(1:B),]
  beta_mode <- apply(b, 2, mode)
  ci <- t(apply(b, 2, quantile, c(0.025, 0.975) ))
  y <- 1:p
  plot(x=beta, y=1:p, ylab="Variable", xlab=expression(beta),
      xlim=range(beta, ci), pch="X", ylim=rev(range(y)), main="Lasso with gamma hyperprior") #LSE with X
  for(i in 1:p) {
    lines(x=ci[i,], y=c(i,i))
    points(x=ci[i,], y=c(i,i), pch="|")
  }
  points(x=beta_mode, y=1:p, col="red") #Lasso posterior mode
  abline(v=0, lty=2)
}

mode <- function(x) {
  dx <- density(x)
  dx$x[which.max(dx$y)]
}
