\name{coef.garchx}
\alias{coef.garchx}
\alias{fitted.garchx}
\alias{logLik.garchx}
\alias{print.garchx}
\alias{residuals.garchx}
\alias{vcov.garchx}

\title{Extraction functions for 'garchx' objects}

\description{Extraction functions for objects of class 'garchx'}

\usage{
  \method{coef}{garchx}(object, ...)
  \method{fitted}{garchx}(object, ...)
  \method{logLik}{garchx}(object, ...)
  \method{print}{garchx}(x, ...)
  \method{residuals}{garchx}(object, ...)
  \method{vcov}{garchx}(object, ...)
}

\arguments{
  \item{object}{an object of class 'garchx'}
  \item{x}{an object of class 'garchx'}
  \item{\dots}{additional arguments}
}

\value{
  \item{coef:}{a numeric vector containing parameter estimates}
  \item{fitted:}{a \code{\link{zoo}} object with the fitted conditional variance}
  \item{logLik:}{log-likelihood (normal density)}
  \item{print:}{a print of the estimation results}
  \item{residuals:}{a \code{\link{zoo}} object with the standardised residuals}
  \item{vcov:}{coefficient variance-covariance matrix}
}

\author{
  Genaro Sucarrat, \url{http://www.sucarrat.net/}
}

\seealso{\code{\link{garchx}}, \code{\link{garchxSim}}, \code{\link{zoo}}}

\examples{
##simulate from a garch(1,1):
set.seed(123)
y <- garchxSim(1000)

##estimate garch(1,1) model:
mymod <- garchx(y)
mymod

##extract coefficients:
coef(mymod)

##extract conditional variances and plot them:
sigma2hat <- fitted(mymod)
plot(sigma2hat)

##extract log-likelihood:
logLik(mymod)

##print results:
print(mymod)

##extract standardised residuals and plot them:
etahat <- residuals(mymod)
plot(etahat)

##extract coefficient variance-covariance matrix:
vcov(mymod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
