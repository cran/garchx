garchx <-
function(y, order=c(1,1), arch=NULL, garch=NULL, asym=NULL,
  xreg=NULL, initial.values=NULL, backcast.values=NULL,
  lower=0, upper=+Inf, control=list(), hessian.control=list(),
  solve.tol=.Machine$double.eps, c.code=FALSE, penalty.value=NULL,
  sigma2.min=.Machine$double.eps, objective.fun=1, verbose=TRUE)
{
  ##create auxiliary list, date, parnames:
  aux <- list()
  aux$date <- date()
  parnames <- "intercept"

  ##y related:
  aux$y.name <- deparse(substitute(y))
  y <- na.trim(as.zoo(y))
  aux$y.n <- NROW(y)
  aux$y.coredata <- coredata(y)
  aux$y.index <- index(y)
  aux$y2 <- aux$y.coredata^2
  aux$y2mean <- mean(aux$y2)

  ##order argument
  ##--------------

  aux$order <- c(0,0,0)
  if( length(order)>0 ){
    aux$order[ 1:length(order) ] <- order[ 1:length(order) ]
  }
  if(is.null(arch) && aux$order[1]>0){ arch <- 1:aux$order[1] }
  if(is.null(garch) && aux$order[2]>0){ garch <- 1:aux$order[2] }
  if(is.null(asym) && aux$order[3]>0){ asym <- 1:aux$order[3] }
  
  ##arch, garch, asym arguments
  ##---------------------------
  
  ##arch:
  if(is.null(arch)){
    aux$archK <- 0
  }else{
    if( identical(arch,0) ){ arch <- NULL }  
    aux$archK <- length(arch)
  }
  aux$arch <- arch

  ##garch:
  if(is.null(garch)){
    aux$garchK <- 0
  }else{
    if( identical(garch,0) ){ garch <- NULL }  
    aux$garchK <- length(garch)
  }
  aux$garch <- garch
  aux$garchOrder <- ifelse(is.null(aux$garch), 0, max(aux$garch))
    
  ##asym:
  if(is.null(asym)){
    aux$asymK <- 0
  }else{
    if( identical(asym,0) ){ asym <- NULL }  
    aux$asymK <- length(asym)
  }
  aux$asym <- asym

  ##xregK, maxpqr, maxpqrpluss1:
  aux$xregK <- ifelse(is.null(xreg), 0, NCOL(xreg))
  aux$maxpqr <- max(aux$archK,aux$garchK,aux$asymK)
  aux$maxpqrpluss1 <- aux$maxpqr + 1
  
  ##parameter indices and names
  ##---------------------------
  
  ##arch:
  if( aux$archK>0 ){
    aux$archIndx <- 2:c( length(arch)+1)
    parnames <- c(parnames, paste0("arch", aux$arch))
  }else{ aux$archIndx <- 0 }
  
  ##garch:
  if( aux$garchK>0 ){
    aux$garchIndx <-
      c( max(1,aux$archIndx) +1) :
      c( max(1,aux$archIndx)+ length(aux$garch) )
    parnames <- c(parnames, paste0("garch", aux$garch))
  }else{ aux$garchIndx <- 0 }
  
  ##asym:
  if( aux$asymK>0 ){
    aux$asymIndx <-
      c( max(1,aux$archIndx,aux$garchIndx) +1) : 
      c( max(1,aux$archIndx,aux$garchIndx) + length(aux$asym) )
    parnames <- c(parnames, paste0("asym", aux$asym))
  }else{ aux$asymIndx <- 0 }
  
  ##xreg:
  if( aux$xregK>0 ){
    aux$xregIndx <-
      c( max(1,aux$archIndx,aux$garchIndx,aux$asymIndx) +1) : 
      c( max(1,aux$archIndx,aux$garchIndx,aux$asymIndx) +aux$xregK)
    xregNames <- colnames(xreg)
    if(is.null(xregNames)){
      xregNames <- paste0("xreg", 1:aux$xregK)
    }
    alreadyTaken <- c(1:aux$xregK)[ xregNames %in% parnames ]
    if( length(alreadyTaken)>0 ){
      xregNames[ alreadyTaken ] <- ""
    }
    missingNames <- which( xregNames=="" )
    if( length(missingNames) > 0 ){
      for(i in missingNames){
        xregNames[i] <- paste0("xreg", i)
      }
    }
    parnames <- c(parnames, xregNames)
  }else{ aux$xregIndx <- 0 }

  
  ##auxiliary vectors and matrices
  ##------------------------------

  ##short y2, short ynotzero
  aux$y2short <- aux$y2[ aux$maxpqrpluss1:aux$y.n ]
  if(objective.fun==0){
    aux$ynotzero <-
      as.numeric(aux$y.coredata != 0)[ aux$maxpqrpluss1:aux$y.n ]
  }

  ##arch matrix:
  if( aux$archK>0 ){
    backvals <- ifelse(is.null(backcast.values),
      aux$y2mean, backcast.values)
    aux$y2matrix <- matrix(backvals, aux$y.n, aux$archK)
    for(i in 1:aux$archK){
      aux$y2matrix[ c(1+aux$arch[i]):aux$y.n ,i] <-
        aux$y2[ 1:c(aux$y.n-aux$arch[i]) ] 
    }
  }

  ##garch vector:
  if( aux$garchK>0 ){
    backvals <- ifelse(is.null(backcast.values),
      aux$y2mean, backcast.values)
    aux$sigma2 <- rep(backvals, aux$y.n)
  }

  ##asym matrix:    
  if( aux$asymK>0 ){
    if(is.null(backcast.values)){
      aux$Ineg <- as.numeric(aux$y.coredata < 0)
      aux$Inegy2 <- aux$Ineg*aux$y2
      backvals <- mean(aux$Inegy2)
    }else{ backvals <- backcast.values }
    aux$Inegy2matrix <- matrix(backvals, aux$y.n, aux$asymK)
    for(i in 1:aux$asymK){
      aux$Inegy2matrix[ c(1+aux$asym[i]):aux$y.n ,i] <-
        aux$Inegy2[ 1:c(aux$y.n-aux$asym[i]) ] 
    }
  }

  ##xreg:  
  if(aux$xregK>0){
    xreg <- na.trim(as.zoo(xreg))
    xreg <-
      window(xreg, start=aux$y.index[1], end=aux$y.index[aux$y.n])
    xreg <- coredata(xreg)
    colnames(xreg) <- NULL
    aux$xreg <- cbind(xreg)
  }

 
  ##initial parameter values:
  ##-------------------------
  if(is.null(initial.values)){

    ##intercept:
    aux$initial.values <- 0.1

    ##arch:
    if( aux$archK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.1/aux$archK, aux$archK)/aux$arch)
#        c(aux$initial.values, rep(0.1/aux$archK, aux$archK))
    }

    ##garch:
    if( aux$garchK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.7/aux$garchK, aux$garchK)/aux$garch)
#        c(aux$initial.values, rep(0.7/aux$garchK, aux$garchK))
    }

    ##asym:
    if( aux$asymK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.02, aux$asymK)/aux$asym)
#        c(aux$initial.values, rep(0.02, aux$asymK))
    }

    ##xreg:
    if( aux$xregK>0 ){
      aux$initial.values <-
        c(aux$initial.values, rep(0.01, aux$xregK))
    }
  }
  ##here, I could add a check for stability; if unstable,
  ##stop or modify?


  ##miscellaneous
  ##-------------

  aux$backcast.values <- backcast.values
  aux$upper <- upper
  aux$lower <- lower  
  aux$control <- control
  aux$hessian.control <- hessian.control
  aux$solve.tol <- solve.tol
  aux$c.code <- c.code
  aux$sigma2.min <- sigma2.min
  aux$objective.fun <- objective.fun
#  aux$require.stability <- require.stability
  aux$penalty.value <- penalty.value
  if(is.null(aux$penalty.value)){
    aux$penalty.value <- garchxObjective(aux$initial.values, aux)
  }

  ##estimation
  ##----------
  
  result <- nlminb(aux$initial.values, garchxObjective,
    aux=aux, control=aux$control, upper=aux$upper, lower=aux$lower)
  names(result$par) <- parnames
   
  ##verbose output?
  ##---------------
  if(verbose){                               

    ##fitted values, residuals:
    sigma2 <- garchxRecursion(as.numeric(result$par), aux)
    residStd <- aux$y.coredata/sqrt(sigma2)
    ##convert to zoo:
    sigma2 <- zoo(sigma2, order.by=aux$y.index)
    residStd <- zoo(residStd, order.by=aux$y.index)
    ##shorten vectors, add to result:
    result$fitted <- sigma2[ aux$maxpqrpluss1:aux$y.n ]
    result$residuals <- residStd[ aux$maxpqrpluss1:aux$y.n ]    

    ##hessian:    
    result$hessian <- optimHess(result$par, garchxObjective,
      aux=aux, control=aux$hessian.control)  

    ##kappahat:
    if(aux$objective.fun==0){
      kappahat <- sum( aux$ynotzero*result$residuals^4)/sum(aux$ynotzero)
    }
    if(aux$objective.fun==1){
      kappahat <- mean( result$residuals^4 )
    }

    ##vcov: divide by nobs for finite sample version
    nobs <- length(result$residuals)
    result$vcov <-
      (kappahat - 1)*solve(result$hessian, tol=aux$solve.tol)/nobs

  }

  ##return result:
  result <- c(aux, result)
  class(result) <- "garchx"
  return(result)

}
