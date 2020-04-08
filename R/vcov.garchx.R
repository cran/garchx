vcov.garchx <-
function(object, ...){
  if(is.null(object$vcov)){
    ##hessian:    
    object$hessian <- optimHess(object$par, garchxObjective,
      object=object, control=object$hessian.control)  
    ##kappahat:
    object$residuals <- residuals.garchx(object)
    if(object$objective.fun==0){
      kappahat <- sum( object$ynotzero * object$residuals^4)/sum(object$ynotzero)
    }
    if(object$objective.fun==1){
      kappahat <- mean( object$residuals^4 )
    }
    ##vcov: divide by n for finite sample version
    nobs <- length(object$residuals)
    object$vcov <-
      (kappahat - 1)*solve(object$hessian, tol=object$solve.tol)/nobs
  }
  return(object$vcov) 
}
