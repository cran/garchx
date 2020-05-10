garchxSim <-
function(n, intercept=0.2, arch=0.1, garch=0.8,
  asym=NULL, xreg=NULL, innovations=NULL, backcast.values=list(),
  verbose=FALSE, as.zoo=TRUE)
{
  ##orders:
  archOrder <- length(arch)
  garchOrder <- length(garch)
  asymOrder <- length(asym)
  maxOrder <- max(archOrder, garchOrder, asymOrder)  
  
  ##initiate:
  if(is.null(innovations)){ innovations <- rnorm(n) }
  z2 <- innovations^2
  Ineg <- as.numeric(innovations < 0)
  sigma2 <- rep(0,n)
  if(is.null(xreg)){ xreg <- rep(0,n) }
  
  ##backcast values
  ##===============
  if(maxOrder > 0){

    ##innovations, Ineg, z2:
    if(is.null(backcast.values$innovations)){
      backcast.values$innovations <- rep(0, maxOrder) 
    }
    innovations <- c(backcast.values$innovations, innovations)
    
    ##z2:
    if(is.null(backcast.values$z2)){
      z2mean <- mean(z2)
      backcast.values$z2 <- rep(z2mean, maxOrder) 
    }
    z2 <- c(backcast.values$z2, z2)

    ##Ineg:
    if(is.null(backcast.values$Ineg)){
      backcast.values$Ineg <- rep(0, maxOrder) 
    }
    Ineg <- c(backcast.values$Ineg, Ineg)

    ##sigma2:
    if(is.null(backcast.values$sigma2)){
      Esigma2 <- intercept/(1-sum(arch)-sum(garch))
      if( abs(Esigma2)==Inf ){ stop("Initial values of sigma2 are not finite") }
      backcast.values$sigma2 <- rep(Esigma2, maxOrder) 
    }
    sigma2 <- c(backcast.values$sigma2, sigma2)

    ##xreg:
    if(is.null(backcast.values$xreg)){
      xregmean <- mean(xreg)
      backcast.values$xreg <- rep(xregmean, maxOrder) 
    }
    xreg <- c(backcast.values$xreg, xreg)
    
  }

  ##recursion
  ##=========
  
  xregsum <- intercept + xreg
  archsum <- garchsum <- asymsum <- 0
  for(i in c(1+maxOrder):length(sigma2) ){

    if(archOrder > 0){
      archsum <-
        sum(arch*z2[c(i-1):c(i-archOrder)]*sigma2[c(i-1):c(i-archOrder)])
    }

    if(garchOrder > 0){
      garchsum <- sum(garch*sigma2[c(i-1):c(i-garchOrder)])
    }

    if(asymOrder > 0){
      asymsum <- sum(asym*Ineg[c(i-1):c(i-asymOrder)]*z2[c(i-1):c(i-asymOrder)]*sigma2[c(i-1):c(i-asymOrder)])
    }

    sigma2[i] <- archsum + garchsum + asymsum + xregsum[i]

  }


  ##prepare result
  ##==============
  
  if(verbose){
    sigma <- sqrt(sigma2)
    y <- sigma*innovations
    result <- cbind(y, sigma, sigma2, Ineg, innovations)
    if(maxOrder > 0){ result <- result[-c(1:maxOrder),] }
  }else{
    sigma <- sqrt(sigma2)
    result <- sigma*innovations
    if(maxOrder > 0){ result <- result[-c(1:maxOrder)] }
  }

  ##return result
  ##=============
  
  if(as.zoo){ result <- as.zoo(result) }
  return(result)

}
