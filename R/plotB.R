#' @title Function to plot a graphical summary of regularized beta coefficients in linear regression
#'
#' @description Function generates plots of least squares estimate, maximum-a-posteriori estimate and 95 percent confidence intervals
#' of beta coefficients for each predictor. Maximum-a-posteriori estimate and 95 percent confidence intervals are generated using
#' the output from either bayesLassoGamma or bayesHorseshoe function.
#'
#' @param x Matrix of predictors, size n by p. Should be centered and scaled before calling the function.
#'
#' @param y Matrix of response, size p by 1. Should be centered before calling the function.
#'
#' @param Total Total number of iterations for Gibbs sampler
#'
#' @param B Burn-in iterations for Gibbs sampler
#'
#' @param thin Thinning parameter, only used if bayesHorseshoe function is specified as an argument
#'
#' @param title Title of the produced plot
#'
#' @param color Color with which maximum-a-posteriori estimates are labeled on the plot
#'
#' @param includeCI Can be TRUE or FALSE, indicates whether to plot 95 percent confidence intervals
#'
#' @param includeMAP Can be TRUE or FALSE, indicates whether to plot maximum-a-posteriori estimates
#'
#' @param includeLSE Can be TRUE or FALSE, indicates whether to plot least squares estimates
#'
#' @param method Can be either bayesLassoGamma or bayesHorseshoe, depending on which regularization method is desired.
#' bayesGammaLasso performs bayesian lasso selection method with gamma hyperprior on \eqn{\lambda^2} as described
#' in Park and Casella (2008). bayesHorseshoe performs horseshoe regularization method as described in Makalic and Schmidt (2016).
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
#' @references
#' Park T, & Casella G. (2008). The Bayesian Lasso. Journal of the American Statistical Association,
#' 103:482, 681-686, DOI: 10.1198/016214508000000337
#'
#' E. Makalic and D. F. Schmidt. (2016) A Simple Sampler for the Horseshoe Estimator.
#' IEEE Signal Processing Letters, 23(1), 179-182. Jan. 2016, doi: 10.1109/LSP.2015.2503725
#'
#' @export plotB

plotB <- function(x,y,Total=10000,B=5000, thin=1,title="Plot", color="red", includeCI=TRUE, includeMAP=TRUE,
                  includeLSE=TRUE, method=c("bayesLassoGamma","bayesHorseshoe")){
  #find LSE
  xtx <- t(x) %*% x
  xy <- t(x) %*% y
  beta <- solve(xtx, xy)

  # find MAP
  if (method=="bayesLassoGamma"){
    b <- bayesLassoGamma(x,y,Total,B)
  }
  if (method=="bayesHorseshoe"){
    b <- bayesHorseshoe(x,y,Total,B,thin)
  }
  p <- ncol(x)
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
