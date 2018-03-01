.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Copyright (c) 2017-2018 Metin Bulus \nSome rights reserved.")
  if(exists(c("plot.mdes", "plot.power"))) {
    packageStartupMessage("------------------------------------ \nSome plotting functions might have been masked by 'PowerUpR'. \nAccess plotting functions directly as 'cosa::plot()'")
  }
}
