\name{garchx}
\alias{garchx}

\title{Estimate a GARCH-X model}

\description{
Normal Maximum Likelihood (ML) estimation of a GARCH(p,q,r)-X model. The estimation procedure will, in general, provide consistent estimates when the data are not normally distributed. A wide range of specifications can be estimated. 
}

\usage{
garchx(y, order = c(1, 1), arch = NULL, garch = NULL, asym = NULL,
  xreg = NULL, initial.values = NULL, backcast.values = NULL,
  lower = 0, upper = +Inf, control = list(), hessian.control = list(),
  solve.tol = .Machine$double.eps, c.code = FALSE, penalty.value = NULL,
  sigma2.min = .Machine$double.eps, objective.fun = 1, verbose = TRUE)
}

\arguments{
  \item{y}{numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{order}{numeric vector of length 1, 2 or 3, for example \code{c(1,1,1)}. The first entry controls the ARCH order, the second the GARCH order and the third the ASYM (asymmetry/leverage) order}
  \item{arch}{\code{NULL} or numeric vector containing the ARCH-terms to include. Note: If not \code{NULL}, then the value of the \code{ARCH} argument overrides the value of the first entry in the \code{order} argument}
  \item{garch}{\code{NULL} or numeric vector containing the GARCH-terms to include. Note: If not \code{NULL}, then the value of the \code{GARCH} argument overrides the value of the second entry in the \code{order} argument}
  \item{asym}{\code{NULL} or numeric vector containing the ASYM-terms (asymmetry/leverage terms) to include. Note: If not \code{NULL}, then the value of the \code{ASYM} argument overrides the value of the third entry in the \code{order} argument}
  \item{xreg}{numeric vector, time-series or \code{\link{zoo}} object. Missing values in the beginning and at the end of the series is allowed, as they are removed with the \code{\link{na.trim}} command}
  \item{initial.values}{\code{NULL} or a numeric vector with the initial parameter values passed on to the optimisation routine, \code{\link{nlminb}}. If \code{NULL}, the default, then the values are chosen automatically}
  \item{backcast.values}{\code{NULL} or 0. The backcast values are used to initiate the forward recursion. If \code{NULL} (default), then the values are chosen automatically. If 0, then the initial recursion values are all set to 0}
  \item{lower}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}}}
  \item{upper}{numeric vector, either of length 1 or the number of parameters to be estimated, see \code{\link{nlminb}}}
  \item{control}{a \code{\link{list}} passed on to the \code{control} argument of \code{\link{nlminb}}}
  \item{hessian.control}{a \code{\link{list}} passed on to the \code{control} argument of \code{\link{optimHess}}}
  \item{solve.tol}{numeric value passed on to the \code{tol} argument of \code{\link{solve}}, which is called whenever the coefficient variance-coariance matrix is computed. The value controls the toleranse for detecting linear dependence between columns when inverting a matrix}
  \item{c.code}{logical, if \code{TRUE} then compiled C code is used in the forward recursion (not available yet)}
  \item{penalty.value}{\code{NULL} (default) or a numeric value. If \code{NULL}, then the log-likelihood value associated with the initial values is used. Sometimes estimation can result in \code{NA} and/or \code{+/-Inf} values. The \code{penalty.value} is the value returned by the objective function \code{\link{garchxObjective}} in the presence of \code{NA} or \code{+/-Inf} values}
  \item{sigma2.min}{numeric}
  \item{objective.fun}{numeric, either 1 or 0}
  \item{verbose}{logical. If code{TRUE} (default), then the coefficient variance-covariance is computed during estimation, and the fitted values and residuals are attached to the returned object. If \code{FALSE}, then these operations are skipped, and hence estimation is faster}
}

\value{A list of class 'garchx'}

\references{
Christian Francq and Le Quien Thieu (2018): 'QML inference for volatility models with covariates', Econometric Theory, doi:10.1017/S0266466617000512
Christian Francq and Jean-Michel Zakoian (2019): 'GARCH Models', 2nd Edition, Wiley
}

\author{Genaro Sucarrat, \url{http://www.sucarrat.net/}}

\seealso{
\code{\link{garchxSim}}, \code{\link{nlminb}}, \code{\link{optimHess}}
}

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

##extract variance-covariance matrix of mean equation:
vcov(mymod)

}

\keyword{Statistical Models}
\keyword{Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
