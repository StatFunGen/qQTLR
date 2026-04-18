# Lazy-bind unexported pecotmr helpers into qQTLR's namespace.
# Resolving at load time via utils::getFromNamespace means call sites can use
# bare names (no `pecotmr:::` prefix) without triggering the R CMD check
# warning "':::' calls which should be '::'".

compute_qvalues        <- NULL
pval_cauchy            <- NULL
drop_collinear_columns <- NULL

.onLoad <- function(libname, pkgname) {
  ns_self <- asNamespace(pkgname)
  for (sym in c("compute_qvalues", "pval_cauchy", "drop_collinear_columns")) {
    assign(sym, utils::getFromNamespace(sym, "pecotmr"), envir = ns_self)
  }
}
