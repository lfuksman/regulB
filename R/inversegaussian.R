#' @title Function to sample from inverse gaussian distribution
#'
#' @description Function intakes a vector of means and a single shape parameter.
#' It then generates a vector of random variables from inverse gaussian distribution
#' from respective mean and same shape parameters.
#'
#' @param mean A vector of mean parameters
#'
#' @param shape Single parameter
#'
#' @return A vector of random variables from inverse gaussian distribution
#' @examples
#'inversegaussian(60,30)
#'inversegaussian(c(35,20),5)
#' @export inversegaussian
#' @references  Michael, John R.; Schucany, William R.; Haas, Roy W. (1976),
#' "Generating Random Variates Using Transformations with Multiple Roots",
#' The American Statistician, 30 (2): 88–90, doi:10.1080/00031305.1976.10479147, JSTOR 2683801

inversegaussian <- function(mean, shape){
  p <- length(mean)
  rinv <- numeric(p)
  v <- rnorm(p,0,1)
  y <- v^2
  x <- mean+(mean^2*y)/(2*shape)-(mean/(2*shape))*sqrt(4*mean*shape*y+mean^2*y^2)
  z <- runif(p,0,1)
  for (i in 1:p){
    if (z[i] <= mean[i]/(mean[i]+x[i])) {
      rinv[i] <- x[i]
    } else{
      rinv[i] <- mean[i]^2/x[i]
    }
  }
  return(rinv)
}
