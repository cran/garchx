quantile.garchx <-
function(x, probs=0.025, names=TRUE,
  type=7, as.zoo=TRUE, ...)
{
  etahat <- residuals.garchx(x)
  sigma <- sqrt(fitted.garchx(x)) 
  qvals <- quantile(etahat, probs=probs, names=names, type=type) 
  iN <- NROW(etahat)
  iCols <- length(qvals)
  result <- matrix(NA, iN, iCols)
  colnames(result) <- names(qvals)
  for(i in 1:iCols){
    result[,i] <- sigma*qvals[i]
  }
  if(iCols==1){ result <- as.vector(result) }
  if(as.zoo){ result <- zoo(result, order.by=index(etahat)) }
  return(result)
}
