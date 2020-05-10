.onAttach <- function(libname, pkgname)
{
  txt <- c("\n",
    paste(sQuote("garchx"), "version 1.1\n"),
    "\n",
    paste0("Flexible and Robust GARCH-X modelling"),
    "\n",
    paste("CRAN website: https://CRAN.R-project.org/package=garchx"),
    "\n")
  
  ##print message:
  if(interactive() || getOption("verbose")){
    packageStartupMessage(paste(strwrap(txt, indent = 2,
      exdent = 4), collapse = "\n"))
  }
} #close .onAttach
