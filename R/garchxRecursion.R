garchxRecursion <-
function(pars, aux)
{
  ##initiate/intercept:
  innov <- rep(pars[1], aux$y.n)
  
  ##arch:
  if( aux$archK>0 ){
    ##use crossprod() instead of %*%?
    innov <- innov + aux$y2matrix %*% pars[ aux$archIndx ]
  }

  ##asym:
  if( aux$asymK>0 ){
    ##use crossprod() instead of %*%?
    innov <- innov + aux$Inegy2matrix %*% pars[ aux$asymIndx ]
  }

  ##xreg:
  if( aux$xregK>0 ){
    ##use crossprod() instead of %*%?
    innov <- innov + aux$xreg %*% pars[ aux$xregIndx ]
  }
    

  ##garchK=0
  ##--------
  if( aux$garchK == 0 ){
    sigma2 <- as.vector(innov)
  }


  ##garchK>0
  ##--------  
  if( aux$garchK > 0 ){
  
    innov <- as.vector(innov)
    sigma2 <- aux$sigma2
    parsgarch <- rep(0, aux$garchOrder)
    parsgarch[ aux$garch ] <- pars[ aux$garchIndx ]
        
    ##if(c.code):
    if(aux$c.code){

#      ##development version:
#      tmp <- GARCHXRECURSION(as.integer(aux$garchOrder),
#        as.integer(aux$recursion.n), as.integer(aux$garchOrder),
#        as.numeric(sigma2), as.numeric(parsgarch),
#        as.numeric(innov))

      ##package version:
      tmp <- .C("GARCHXRECURSION", iStart = as.integer(aux$garchOrder),
        iEnd = as.integer(aux$recursion.n), iGARCHorder = as.integer(aux$garchOrder),
        sigma2 = as.double(sigma2), parsgarch = as.double(parsgarch),
        innov = as.double(innov), PACKAGE = "garchx")

      ##sigma2:
      sigma2 <- tmp$sigma2

    }else{

      ##if(garch1):
      if(aux$garchOrder==1){
        for( i in 2:aux$y.n ){
          sigma2[i] <- parsgarch * sigma2[i-1] + innov[i]
        }
      }else{
        for( i in c(1+aux$garchOrder):aux$recursion.n ){
#OLD:
#        for( i in c(1+aux$garchOrder):aux$y.n ){
          sigma2[i] <-
            sum(parsgarch*sigma2[ c(i-1):c(i-aux$garchOrder) ]) + innov[i]
        }
      }

    } #close if(c.code)

  } #close if(garchK==0)else(..)
  

  ##output:
  return(sigma2)
  ##note: this volatility contains the first observations
  ##(i.e. 1:aux$maxpqr), which should be deleted before
  ##returnted by fitted(..) etc.

}
