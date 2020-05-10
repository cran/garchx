garchxObjective <-
function(pars, aux)
{
  ##initiate:
  parsOK <- TRUE; sigma2OK <- TRUE
  
  ##check parameters:
  if( any(is.na(pars))){ parsOK <- FALSE  }
#  if( any(is.na(pars)) || any(pars<aux$lower) || any(pars>aux$upper) ){
#    parsOK <- FALSE
#  }
#  if( parsOK && aux$require.stability ){
#    parsOK <- sum(pars[2:3]) < 1
#  }

  ##compute and check sigma2:
  if( parsOK ){
    sigma2 <- garchxRecursion(pars, aux)
    sigma2 <- sigma2[ aux$maxpqrpluss1:aux$y.n ]
    if( any(sigma2<=0) || any(sigma2==Inf) ){
      sigma2OK <- FALSE
    }else{
      ##robustify log-transform in log-likelihood:
      sigma2 <- pmax( aux$sigma2.min, sigma2 )
    }
  }

  ##compute average log-likelihood:
  if( parsOK && sigma2OK ){
    if(aux$objective.fun==0){
      result <- mean( aux$ynotzero*(aux$y2short/sigma2 + log(sigma2)) )
    }
    if(aux$objective.fun==1){
      result <- mean( aux$y2short/sigma2 + log(sigma2) )
    }  
  }else{
    result <- aux$penalty.value
  }

  ##return result:
  return(result)

}
