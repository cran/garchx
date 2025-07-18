\name{garchx-package}
\alias{garchx-package}
\docType{package}

\title{Flexible and Robust GARCH-X Modelling}

\description{
Flexible and robust estimation and inference of GARCH(q,p,r)-X models, where q is the GARCH order, p is the ARCH order, r is the asymmetry or leverage order, and 'X' indicates that covariates can be included. Suitable subsets of the coefficients can be restriced to zero by omission, and Quasi Maximum Likelihood (QML) methods ensure estimates are generally consistent, even when the standardised innovations are non-normal and/or dependent.
}

\details{
  \tabular{ll}{
  Package: \tab garchx\cr
  Type: \tab Package\cr
  Version: \tab 1.6\cr
  Date: \tab 2025-07-09\cr
  License: \tab GPL-2\cr
  }
}

\author{
Genaro Sucarrat, \url{https://www.sucarrat.net/}\cr\cr
Maintainer: Genaro Sucarrat
}

\seealso{
\code{\link{garchxSim}}, \code{\link{coef}}, \code{\link{fitted}}, \code{\link{logLik}}, \code{\link{print}}, \code{\link{residuals}}, \code{\link{vcov}}
}

\examples{
##simulate from a garch(1,1):
set.seed(123)
y <- garchxSim(1000)

##estimate garch(1,1) model:
mymod <- garchx(y)
mymod

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
