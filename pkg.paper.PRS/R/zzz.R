################################################################################

.onLoad <- function(libname, pkgname) {
  options(
    nboot = 1e5
  )
}

################################################################################

.onUnload <- function(libpath) {
  options(
    nboot = NULL
  )
}

################################################################################
