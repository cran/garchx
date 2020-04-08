fitted.garchx <-
function(object, ...){
  if(is.null(object$fitted)){
    object$fitted <- garchxRecursion(as.numeric(object$par), object)
    object$fitted <- zoo(object$fitted, order.by=object$y.index)
    object$fitted <- object$fitted[ object$maxpqrpluss1:object$y.n ]
  }
  return(object$fitted)
}
