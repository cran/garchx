logLik.garchx <-
function(object, ...){
  y <- object$y.coredata[ object$maxpqrpluss1:object$y.n ]
  sigma2 <- coredata(fitted.garchx(object))
  dnormvals <- dnorm(y, mean=0, sd=sqrt(sigma2), log=TRUE)
  if(object$objective.fun==0){
    result <- sum( object$ynotzero * dnormvals )
  }
  if(object$objective.fun==1){
    result <- sum(dnormvals)
  }
  attr(result, "df") <- length(object$par)
  attr(result, "nobs") <- ifelse(object$objective.fun==1,
    length(sigma2), sum(object$ynotzero))
  return(result)
}
