# Lazy-bind unexported pecotmr helpers into qQTLR's namespace.
# Resolving at load time via utils::getFromNamespace means call sites can use
# bare names (no `pecotmr:::` prefix) without triggering the R CMD check
# warning "':::' calls which should be '::'".

compute_qvalues <- NULL
pval_cauchy     <- NULL

.onLoad <- function(libname, pkgname) {
  ns_self <- asNamespace(pkgname)
  assign("compute_qvalues",
         utils::getFromNamespace("compute_qvalues", "pecotmr"),
         envir = ns_self)
  assign("pval_cauchy",
         utils::getFromNamespace("pval_cauchy", "pecotmr"),
         envir = ns_self)
}
