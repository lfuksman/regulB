\name{diabetes}
\alias{diabetes}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{Diabetes dataset
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("diabetes")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(diabetes)
## maybe str(diabetes) ; plot(diabetes) ...
}
\keyword{datasets}

save("diabetes.txt", file="data/diabetes.RData")
