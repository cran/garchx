ttest0 <-
function(x, k=NULL)
{
  if(class(x)!="garchx"){ stop("'x' not of class 'garchx'") }
  
  ##prepare:
  coefs <- coef.garchx(x)
  coefsNames <- names(coefs)
  n <- length(coefs)
  if(is.null(k)){ k <- 2:n }
  mSigmahat <- vcov.garchx(x)
  
  ##make matrix w/tests:
  result <- matrix(NA, length(k), 4)
  colnames(result) <- c("coef", "std.error", "t-stat", "p-value")
  rownames(result) <- coefsNames[k]
  result[,"coef"] <- coefs[k]
  
  ##fill matrix w/tests:
  for(i in 1:NROW(result)){
  
    evector <- rep(0, n)
    evector[ k[i] ] <- 1
    stderror <-
      as.numeric(sqrt( rbind(evector) %*% mSigmahat %*% cbind(evector) ))
    statistic <- coefs[k[i]]/stderror
    
    result[i,"std.error"] <- stderror
    result[i,"t-stat"] <- statistic
    result[i,"p-value"] <- pnorm(statistic, lower.tail=FALSE)
  
  }

  ##result:
  return(result)

}
