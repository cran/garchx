garchxAvar <-
function(pars, arch=NULL, garch=NULL, asym=NULL,
  xreg=NULL, vcov.type=c("ordinary", "robust"),  innovations=NULL,
  Eeta4=NULL, n=1000000, objective.fun=1, seed=NULL)
{
  ##arguments
  ##---------

  ##pars:
  if( length(pars)==0 ){ stop("length(pars) cannot be 0") }

  ##vcov.type:
  types <- c("ordinary", "robust")
  whichType <- charmatch(vcov.type[1], types)
  vcov.type <- types[whichType]
  if( vcov.type=="robust"){ stop("Sorry, not implemented yet!") }

  ##initiate
  ##--------
  parnames <- "intercept"
  aux <- list()
  aux$y.n <- n
  aux$recursion.n <- aux$y.n #used in recursion only; needed for robust vcov
  aux$y.coredata <- 1 #preliminary placeholder
  aux$y2 <- aux$y.coredata^2
  aux$y2mean <- 1 #preliminary placeholder

  ##arch, garch, asym arguments
  ##---------------------------

  ## note: the K refers to how many arch/asym/garch terms there are,
  ## NOT the arch/asym/garch order

  ##arch:
  if(is.null(arch)){
    aux$archK <- 0
  }else{
    if( identical(arch,0) ){ arch <- NULL }
    aux$archK <- length(arch)
  }
  aux$arch <- arch
  aux$archOrder <- ifelse(is.null(aux$arch), 0, max(aux$arch))

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
  aux$asymOrder <- ifelse(is.null(aux$asym), 0, max(aux$asym))

  ##xregK, maxpqr, maxpqrpluss1:
  aux$xregK <- ifelse(is.null(xreg), 0, NCOL(xreg))
  aux$maxpqr <- max(aux$archOrder,aux$garchOrder,aux$asymOrder)
  aux$maxpqrpluss1 <- aux$maxpqr + 1

  ##parameter indices and names
  ##---------------------------

  ##arch:
  if( aux$archK>0 ){
    aux$archIndx <- 2:c( length(arch)+1 )
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

  ##simulate
  ##--------

  ##arch argument:
  if( aux$archK==0 ){
    archArg <- NULL
  }else{
    archArg <- rep(0, aux$archOrder)
    archArg[ aux$arch ] <- pars[ aux$archIndx ]
  }

  ##garch argument:
  if( aux$garchK==0 ){
    garchArg <- NULL
  }else{
    garchArg <- rep(0, aux$garchOrder)
    garchArg[ aux$garch ] <- pars[ aux$garchIndx ]
  }

  ##asym argument:
  if( aux$archK==0 ){
    asymArg <- NULL
  }else{
    asymArg <- rep(0, aux$asymOrder)
    asymArg[ aux$asym ] <- pars[ aux$asymIndx ]
  }

  ##xreg argument:
  if( aux$xregK==0 ){
    xregArg <- NULL
  }else{
    xregArg <- cbind(xreg) %*% pars[ aux$xregIndx ]
  }

  ##set simulation length:
  if( aux$garchOrder>0 && is.null(innovations) ){
    nadj <- n+aux$garchOrder
  }else{ nadj <- n }

  ##simulate:
  if( !is.null(seed) ){ set.seed(seed) }
  mY <- garchxSim(nadj, intercept=pars[1], arch=archArg, garch=garchArg,
    asym=asymArg, xreg=xregArg, innovations=innovations, verbose=TRUE)
  if( nadj > n ){
    backcast.values <- as.numeric(mY[1:aux$garchOrder,"sigma2"])
    mY <- mY[ -c(1:aux$garchOrder),]
  }else{
    backcast.values <- NULL
  }

  ##y, y2, y2mean:
  aux$y.coredata <- coredata(mY[,"y"])
  aux$y2 <- aux$y.coredata^2
  aux$y2mean <- mean(aux$y2) #default garch backcast value

  ##auxiliary vectors and matrices
  ##------------------------------

  ##short y2, short ynotzero
  aux$y2short <- aux$y2[ aux$maxpqrpluss1:aux$y.n ]
  if( objective.fun==0 ){
    aux$ynotzero <-
      as.numeric(aux$y.coredata != 0)[ aux$maxpqrpluss1:aux$y.n ]
  }

  ##arch matrix:
  if( aux$archK>0 ){
    aux$y2matrix <- matrix(aux$y2mean, aux$y.n, aux$archK)
    for(i in 1:aux$archK){
      aux$y2matrix[ c(1+aux$arch[i]):aux$y.n ,i] <-
        aux$y2[ 1:c(aux$y.n-aux$arch[i]) ]
    }
  }

  ##garch vector:
  if( aux$garchK>0 ){
    aux$sigma2 <- rep(aux$y2mean, aux$y.n)
    if( !is.null(backcast.values) ){
      aux$sigma2[1:aux$garchOrder] <- backcast.values
    }
  }

  ##asym matrix:
  if( aux$asymK>0 ){
    aux$Ineg <- as.numeric(aux$y.coredata < 0)
    aux$Inegy2 <- aux$Ineg*aux$y2
    backvals <- mean(aux$Inegy2)
    aux$Inegy2matrix <- matrix(backvals, aux$y.n, aux$asymK)
    for(i in 1:aux$asymK){
      aux$Inegy2matrix[ c(1+aux$asym[i]):aux$y.n ,i] <-
        aux$Inegy2[ 1:c(aux$y.n-aux$asym[i]) ]
    }
  }

  ##xreg:
  if(aux$xregK>0){
    colnames(xreg) <- NULL
    aux$xreg <- cbind(xreg)
  }

  ##miscellaneous
  ##-------------

  aux$upper <- Inf
  aux$lower <- 0
  aux$control <- list()
  aux$hessian.control <- list()
  aux$solve.tol <- .Machine$double.eps
  aux$c.code <- TRUE
  aux$sigma2.min <- .Machine$double.eps
  aux$objective.fun <- objective.fun
  aux$penalty.value <- garchxObjective(pars, aux)

  ##hessian
  ##-------

  names(pars) <- parnames
  aux$hessian <- optimHess(pars, garchxObjective, aux=aux,
    control=aux$hessian.control)

  ##ordinary vcov
  ##-------------
  if( vcov.type=="ordinary" ){

    ##kappa:
    if( is.null(Eeta4) ){
      if( aux$objective.fun==1 ){
        Eeta4 <- mean( mY[,"innovations"]^4 )
      }
      if( aux$objective.fun==0 ){
        etanotzero <- as.numeric( mY[,"innovations"] != 0 )
        Eeta4 <- sum( etanotzero*innovations^4 )/sum(etanotzero)
      }
    }

    ##avar:
    result <- (Eeta4 - 1)*solve(aux$hessian, tol=aux$solve.tol)
  }


  ##robust vcov
  ##-----------
  if( vcov.type=="robust" ){
    ##tba
  }

  ##return
  ##------
  return(result)

}
