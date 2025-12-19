.onLoad <- function(libname, pkgname) {
  # Register C-callables (R_RegisterCCallable) for cross-package C++ use.
  # This keeps the user-facing API unchanged while enabling R_GetCCallable
  # access from other packages (e.g. EMC2).
  .Call("_AccumulatR_register_ccallables_cpp", PACKAGE = pkgname)
}

