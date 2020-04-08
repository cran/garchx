.onAttach <- function(libname, pkgname)
{
  ##set start-up message:
  txt <- c("\n",
    paste(sQuote("garchx"), "version 1.0\n"),
    "\n",
    paste0("Flexible and Robust GARCH-X Modelling"),
#    "\n",
#    paste("CRAN website: https://CRAN.R-project.org/package=garchx"),
    "\n")
  
  ##print message:
  if(interactive() || getOption("verbose")){
    packageStartupMessage(paste(strwrap(txt, indent = 2,
      exdent = 4), collapse = "\n"))
  }
} #close .onAttach
