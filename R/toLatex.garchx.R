toLatex.garchx <-
function(object, digits=4, ...)
{
  ##equation:
  ##---------

  coefs <- coef.garchx(object)
  coefsNames <- names(coefs)
  coefsNames[1] <- "" #intercept
  coefs <- as.numeric(coefs)
  stderrs <- as.numeric(sqrt(diag(vcov(object))))

  eqtxt <- NULL
  for(i in 1:length(coefs) ){
    ifpluss <- ifelse(i==1, "", " + ")
    eqtxt <- paste(eqtxt,
      ifelse(coefs[i] < 0, " - ", ifpluss), "\\underset{(",
      format(round(stderrs[i], digits=digits), nsmall=digits), ")}{",
      format(round(abs(coefs[i]), digits=digits), nsmall=digits), "}",
      coefsNames[i], sep="")
  }

  txtAddEq <- " \\\\[1mm]"
  eqtxt <- paste0("  \\widehat{\\sigma}_t^2 &=& ", eqtxt, "",
    txtAddEq, " \n")

  ##goodness of fit:
  ##----------------

  goftxt <- NULL
  goftxt <- "   &&"
  iT <- nobs.garchx(object)
  goftxt <- paste(goftxt, " \\text{Log-likelihood: }",
    format(round(as.numeric(logLik.garchx(object)), digits=digits), nsmall=digits),
      "\\qquad T = ", iT, " \\nonumber \n", sep="")

  ##print code:
  ##-----------

  cat("\\begin{eqnarray}\n")
  cat(eqtxt)
  cat(goftxt)
  cat("\\end{eqnarray}\n")

}
