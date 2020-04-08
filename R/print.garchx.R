print.garchx <-
function(x, ...){
  
  ##out1:
  pars <- coef.garchx(x)
  vcovmat <- vcov.garchx(x)
  out1 <- rbind(pars, sqrt(diag(vcovmat)))
  rownames(out1) <- c("Estimate:", "Std. Error:")

  ##out2:
  out2 <- as.data.frame(matrix(NA, 1, 1))
  out2[1, 1] <- as.character(round(logLik.garchx(x), digits = 4))
  rownames(out2) <- "Log-likelihood:"
  colnames(out2) <- ""
    
  ##print message:
  cat("\n")
  cat("Date:", x$date, "\n")
  cat("Method: Gaussian ML\n")
  cat("Message (nlminb):", x$message, "\n")
  cat("No. of observations:", x$y.n - x$maxpqr, "\n")
  cat("Sample:", as.character(x$y.index[1]), "to", as.character(x$y.index[x$y.n]), 
      "\n")
  cat("\n")
#  cat("Coefficients:\n")
  print(out1)
  print(out2)
  cat("\n")
}
