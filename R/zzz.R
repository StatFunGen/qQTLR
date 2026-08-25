# Lazy-bind unexported or renamed pecotmr helpers into qQTLR's namespace.

compute_qvalues              <- NULL
pval_cauchy                  <- NULL
drop_collinear_columns       <- NULL
build_twas_score_row         <- NULL
parse_variant_id             <- NULL
filter_variants_by_ld_reference <- NULL
ld_prune_by_correlation      <- NULL
ld_clump_by_score            <- NULL

.onLoad <- function(libname, pkgname) {
  ns_self <- asNamespace(pkgname)

  # Helper: bind a renamed pecotmr function under a legacy snake_case name.
  bind_from_pecotmr <- function(new_name, old_name, fallback_msg) {
    fn <- tryCatch(
      utils::getFromNamespace(new_name, "pecotmr"),
      error = function(e) NULL
    )
    if (!is.null(fn)) {
      assign(old_name, fn, envir = ns_self)
    } else {
      msg <- fallback_msg
      assign(old_name, function(...) stop(msg), envir = ns_self)
    }
  }

  # pval_cauchy: Cauchy combination test, renamed to pvalAcat in current pecotmr.
  bind_from_pecotmr("pvalAcat", "pval_cauchy",
    "'pvalAcat' is not available in the installed version of pecotmr")

  # drop_collinear_columns: renamed to dropCollinearColumns in current pecotmr.
  bind_from_pecotmr("dropCollinearColumns", "drop_collinear_columns",
    "'dropCollinearColumns' is not available in the installed version of pecotmr")

  # filter_variants_by_ld_reference: renamed to filterVariantsByLdReference.
  bind_from_pecotmr("filterVariantsByLdReference", "filter_variants_by_ld_reference",
    "'filterVariantsByLdReference' is not available in the installed version of pecotmr")

  # ld_prune_by_correlation: renamed to ldPruneByCorrelation.
  bind_from_pecotmr("ldPruneByCorrelation", "ld_prune_by_correlation",
    "'ldPruneByCorrelation' is not available in the installed version of pecotmr")

  # ld_clump_by_score: renamed to ldClumpByScore.
  bind_from_pecotmr("ldClumpByScore", "ld_clump_by_score",
    "'ldClumpByScore' is not available in the installed version of pecotmr")

  # compute_qvalues: q-value estimation via qvalue package, with FDR fallback.
  assign("compute_qvalues", function(pvals) {
    tryCatch(
      qvalue::qvalue(pvals)$qvalues,
      error = function(e) p.adjust(pvals, method = "fdr")
    )
  }, envir = ns_self)

  # build_twas_score_row: not in current pecotmr; used only in the TWAS
  # pipeline path which is not called when --no-twas-weight-calculate is set.
  assign("build_twas_score_row", function(...) stop(
    "'build_twas_score_row' is not available in the installed version of pecotmr"
  ), envir = ns_self)

  # parse_variant_id was renamed to parseVariantId in pecotmr; provide a
  # snake_case alias that calls the camelCase version, with column-name
  # normalisation so callers expecting $chrom/$pos/$A2/$A1 still work.
  parse_variant_id_fn <- tryCatch(
    utils::getFromNamespace("parseVariantId", "pecotmr"),
    error = function(e) NULL
  )
  if (!is.null(parse_variant_id_fn)) {
    alias <- function(ids) {
      res <- parse_variant_id_fn(ids)
      if (is.data.frame(res) && !"chrom" %in% names(res)) {
        if ("CHR" %in% names(res)) names(res)[names(res) == "CHR"] <- "chrom"
      }
      res
    }
    assign("parse_variant_id", alias, envir = ns_self)
  } else {
    assign("parse_variant_id", function(...) stop(
      "'parseVariantId' is not available in the installed version of pecotmr"
    ), envir = ns_self)
  }
}
