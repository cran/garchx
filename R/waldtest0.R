waldtest0 <-
function(x, r=0, R=NULL, level=c(0.1,0.05,0.01),
  vcov.type=NULL, quantile.type=7, n=20000)
{
  if(class(x)!="garchx"){ stop("'x' not of class 'garchx'") }
  
  ##coefs:
  coefs <- as.numeric(coef.garchx(x))
  nmin1 <- length(coefs)-1
  
  ##combination matrix R:
  if(is.null(R)){
    R <- matrix(0, nmin1, nmin1)
    diag(R) <- 1
    R <- cbind( rep(0,nmin1), R)
  }else{
    R <- as.matrix(R)
  }
  rankR <- qr(R)$rank
  if( rankR != NROW(R)){ stop("'R' does not have full rank") }

  ##modify r?
  if( NROW(R)>length(r) ){ r <- rep(r[1], NROW(R)) }
  r <- cbind(r)
  
  ##statistic
  ##---------
  
  SigmahatT <- vcov.garchx(x, vcov.type=vcov.type)
  term3 <- R %*% coefs - r
  term1 <- t(term3)
  term2 <- solve(R %*% SigmahatT %*% t(R))
  statistic <- as.numeric(term1 %*% term2 %*% term3)
  
  ##critical values
  ##---------------
  
  Sigmahat <- nobs(x)*SigmahatT
  Zhat <- t(rmnorm(n, vcov=Sigmahat))
  RZhat <- matrix(NA, NROW(R), n)
  normRZhat <- rep(NA, n) #norm of RZhat
  for(i in 1:n){
    RZhat[,i] <- R %*% Zhat[,i]
    normRZhat[i] <- sum(RZhat[,i]^2)
  }
  critical.values <-
    quantile(normRZhat, prob=1-level, type=quantile.type)
  names(critical.values) <- paste0(100*level,"%")
    
  ##result
  ##------
  
  result <- list(statistic=statistic, critical.values=critical.values)
  return(result)

}
