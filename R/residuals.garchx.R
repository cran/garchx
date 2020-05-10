residuals.garchx <-
function(object, as.zoo=TRUE, ...){
  if(is.null(object$residuals)){
    sigma2 <- garchxRecursion(as.numeric(object$par), object)
    object$residuals <- object$y.coredata/sqrt(sigma2)
    if(as.zoo){
      object$residuals <- zoo(object$residuals, order.by=object$y.index)
    }
    object$residuals <-
      object$residuals[ object$maxpqrpluss1:object$y.n ]
  }
  return(object$residuals)
}
