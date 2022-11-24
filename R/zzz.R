
.onAttach <- function(libname, pkgname) {
  
  options(
    "pboptions" = list(
      type = if (interactive()) "timer" else "none",
      char = "-",
      txt.width = 50,
      gui.width = 300,
      style = 1,
      initial = 0,
      title = "R progress bar",
      label = "",
      nout = 40L,
      min_time = 10
    )
  )
  
  invisible(NULL)
}

.onUnload <- function(libpath) library.dynam.unload("SAMtool", libpath)
