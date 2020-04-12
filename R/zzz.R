## Set these options when package is loaded

.onLoad <- function(libname = find.package("sleepsimReval"), pkgname="sleepsimReval") {
  # Hard-code options in the package
  options(
    sleepsimReval = list(
      # Future library plan
      "plan" = "sequential"
    )
  )
}
