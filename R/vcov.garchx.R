vcov.garchx <-
function(object, vcov.type=NULL, ...)
{
  ##compute hessian?
  ##----------------
  if( ! ("hessian" %in% names(object)) ){
    object$hessian <- optimHess(object$par, garchxObjective,
      aux=object, control=object$hessian.control)  
  }

  ##determine vcov.type
  ##-------------------
  vcovTypes <- c("ordinary", "robust")
  if( is.null(vcov.type) ){
    sysCall <- as.list(object$sys.call)
    if( "vcov.type" %in% names(sysCall) ){
      whichOne <- which( "vcov.type" == names(sysCall) )
      vcov.type <- sysCall[[ whichOne ]]
    }else{
      vcov.type <- vcovTypes
    }
  }    
  whichType <- charmatch(vcov.type[1], vcovTypes)
  vcov.type <- vcovTypes[ whichType ]

  ##vcov comment
  ##------------
  vcovComment <-
    ifelse("vcov" %in% names(object), comment(object$vcov), "NULL")

  ##ordinary vcov
  ##-------------
  if( vcov.type=="ordinary" && vcovComment != "ordinary" ){

    ##kappahat:
    if( is.null(object$residuals) ){
      object$residuals <- residuals.garchx(object)
    }
    if( object$objective.fun==0 ){
      kappahat <-
        sum( object$ynotzero * object$residuals^4)/sum(object$ynotzero)
    }
    if( object$objective.fun==1 ){
      kappahat <- mean( object$residuals^4 )
    }

    ##vcov:
    iN <- length(object$residuals) #divide by n for finite sample version
    object$vcov <-
      (kappahat - 1)*solve(object$hessian, tol=object$solve.tol)/iN
    comment(object$vcov) <- "ordinary"

  } #close if(ordinary)
        
  ##robust vcov
  ##------------
  if( vcov.type=="robust" && vcovComment!="robust" ){

    ##inverse of J:
    Jinv <- solve(object$hessian, tol=object$solve.tol)
    
    ##temporary function:
    funtmp <- function(i, pars){
      object$recursion.n <- i
      object$par <- pars
      sigma2 <- garchxRecursion(as.numeric(object$par), object)
      return(sigma2[i])
    }

    ##matrix w/gradients
    pars <- as.numeric(object$par)
    mIadj <- matrix(NA, object$y.n, length(object$par))
#    sigma2Check <- rep(NA, object$y.n)
    for(i in 1:object$y.n){
      mIadj[i,] <-
        attr(numericDeriv(quote(funtmp(i,pars)),"pars"), "gradient")
#      ##code that enables a of sigma2:
#      tmp <- numericDeriv(quote(funtmp(i,pars)),"pars")
#      sigma2Check[i] <- tmp[1]
#      mIadj[i,] <- attr(tmp, "gradient")
    }

    ##Iadj:
    sigma2 <- garchxRecursion(as.numeric(object$par), object)
    etahatadj <- object$y2/(sigma2^2)
    mIadj <- etahatadj * mIadj
    mIadj <- mIadj[ object$maxpqrpluss1:object$y.n, ]
    iN <- length(residuals.garchx(object))    
    mIadj <- crossprod(mIadj)/iN - object$hessian

    ##vcov:
    object$vcov <- (Jinv %*% mIadj %*% Jinv)/iN
    comment(object$vcov) <- "robust"

  } #close if("robust")
    
  ##return
  ##------
  return(object$vcov) 

}
