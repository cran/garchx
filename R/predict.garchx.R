predict.garchx <-
function(object, n.ahead=10, newxreg=NULL,
  newindex=NULL, n.sim=5000, verbose=FALSE, ...)
{
  ##coefficients
  coefs <- as.numeric(coef.garchx(object))
  interceptCoef <- coefs[1]

  archCoef <- NULL
  if(object$archK>0){
    archCoef <- rep(0, object$archOrder)
    archCoef[ object$arch ] <- coefs[ object$archIndx ]
  }

  garchCoef <- NULL
  if(object$garchK>0){
    garchCoef <- rep(0, object$garchOrder)
    garchCoef[ object$garch ] <- coefs[ object$garchIndx ]
  }

  asymCoef <- NULL
  if(object$asymK>0){
    asymCoef <- rep(0, object$asymOrder)
    asymCoef[ object$asym ] <- coefs[ object$asymIndx ]
  }

  if( object$xregK>0 ){ xregCoef <- coefs[ object$xregIndx ] }


  ##backcast values
  backcast.values <- list()
  if( object$maxpqr>0 ){
    backcast.values$innovations <-
      tail(as.numeric(residuals.garchx(object)), n=object$maxpqr)
    backcast.values$z2 <- backcast.values$innovations^2
    backcast.values$Ineg <- as.numeric( backcast.values$innovations<0 )
    backcast.values$sigma2 <-
      tail(as.numeric(fitted.garchx(object)), n=object$maxpqr)
    if(object$xregK>0){
      backcast.values$xreg <-
        tail(as.numeric(object$xreg %*% xregCoef), n=object$maxpqr)
    }
  }

  ##newxreg
  xreg <- NULL
  if( object$xregK>0 ){
    if( (NROW(newxreg)!=n.ahead) ){
      stop("NROW(newxreg) is unequal to n.ahead")
    }
    xreg <- cbind(newxreg) %*% xregCoef
  }

  ##bootstrap the innovations
  etahat <- coredata(residuals.garchx(object))
  draws <- runif(n.ahead * n.sim, min = 0.5 + 1e-09,
    max = length(etahat) + 0.5 - 1e-09)
  draws <- round(draws, digits = 0)
  innovations <- matrix(etahat[draws], n.ahead, n.sim)

  ##make predictions
  predictions <- matrix(NA, n.ahead, n.sim)
  for(i in 1:n.sim){
    predictions[,i] <- garchxSim(n.ahead, intercept=interceptCoef,
      arch=archCoef, garch=garchCoef, asym=asymCoef, xreg=xreg,
      innovations=innovations[,i], backcast.values=backcast.values,
      verbose=TRUE, as.zoo=FALSE)[,"sigma2"]
  }

  ##result
  result <- as.vector(rowMeans(predictions))
  if(verbose){
    result <- cbind(result, predictions)
    colnames(result) <- c("sigma2", 1:NCOL(predictions))
  }
  if(is.null(newindex)){ newindex <- 1:n.ahead }
  result <- zoo(result, order.by=newindex)
  return(result)

}
