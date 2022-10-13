#' @title plotB
#'
#' @description Generates plots of least squares estimate, maximum-a-posteriori estimate and 95 percent confidence intervals
#' of beta coefficients for each predictor. Maximum-a-posteriori estimate and 95 percent confidence intervals are generated using
#' the output from bayesLassoGamma function.
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
#' plotB(x,y,color="blue",title="Lasso with Gamma hyperprior",includeMAP=FALSE,includeCI=TRUE,includeLSE =FALSE) #CI
#' plot(x,y,color="red",title="Lasso with Gamma hyperprior",includeMAP=TRUE,includeCI=TRUE,includeLSE =TRUE) #CI,MAP,LSE
#'
#' @export plotB
#' @importFrom
#' @references

plotB <- function(x,y,Total=10000,B=5000, title="Plot", color="red", includeCI=TRUE, includeMAP=TRUE,includeLSE=TRUE){
  #find LSE
  xtx <- t(x) %*% x
  xy <- t(x) %*% y
  beta <- solve(xtx, xy)

  # find MAP
  b <- bayesLassoGamma(x,y,Total,B)
  p <- ncol(x)
  b <- b[-(1:B),]
  beta_MAP <- apply(b, 2, MAP)
  # find CI

  ci <- t(apply(b, 2, quantile, c(0.025, 0.975) ))
  y <- 1:p


  if (includeLSE==TRUE){ #if include LSE
    plot(x=beta, y=1:p, ylab="Variable", xlab=expression(beta),
        xlim=range(beta, ci), pch="X", ylim=rev(range(y)), main=title) #LSE with X
    abline(v=0, lty=2)
    if (includeCI==TRUE){
      for(i in 1:p) {
        lines(x=ci[i,], y=c(i,i))
        points(x=ci[i,], y=c(i,i), pch="|")
      }
    }
    if (includeMAP==TRUE){
      points(x=beta_MAP, y=1:p, col=color) #Lasso posterior mode
    }
  } else{ # if don't include LSE
    plot(x=NA, type="n",ylab="Variable", xlab=expression(beta),
         xlim=range(beta, ci), pch="X", ylim=rev(range(y)), main=title)
    abline(v=0, lty=2)
    if (includeCI==TRUE){
      for(i in 1:p) {
        lines(x=ci[i,], y=c(i,i))
        points(x=ci[i,], y=c(i,i), pch="|")
      }
    }
    if (includeMAP==TRUE){
      points(x=beta_MAP, y=1:p, col=color) #Lasso posterior mode
    }
  }
}

MAP <- function(x) {
  dx <- density(x)
  dx$x[which.max(dx$y)]
}
