residuals.garchx <-
function(object, ...){
  if(is.null(object$residuals)){
    sigma2 <- garchxRecursion(as.numeric(object$par), object)
    object$residuals <- object$y.coredata/sqrt(sigma2)
    object$residuals <- zoo(object$residuals, order.by=object$y.index)
    object$residuals <-
      object$residuals[ object$maxpqrpluss1:object$y.n ]
  }
  return(object$residuals)
}
