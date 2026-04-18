#' Univariate regression with z/p/q augmentation
#'
#' Thin wrapper around \code{susieR::univariate_regression} that strips rows
#' with NA in \code{y}, then augments the returned list with z-scores,
#' p-values, q-values, and names derived from \code{colnames(X)}.
#'
#' @param X Numeric matrix of genotypes (n x p).
#' @param y Numeric vector of phenotypes (length n).
#' @param Z Optional numeric matrix of covariates (n x k).
#' @param center Logical, whether to center (default: TRUE).
#' @param scale Logical, whether to scale (default: FALSE).
#' @param return_residuals Logical, whether to return residuals (default: FALSE).
#' @return A list with \code{betahat}, \code{sebetahat}, \code{z_scores},
#'   \code{p_values}, \code{q_values}, and optionally \code{residuals}.
#' @noRd
univariate_regression <- function(X, y, Z = NULL, center = TRUE,
                                  scale = FALSE, return_residuals = FALSE) {
  y_na <- which(is.na(y))
  if (length(y_na)) {
    X <- X[-y_na, ]
    y <- y[-y_na]
    if (!is.null(Z)) Z <- Z[-y_na, , drop = FALSE]
  }

  nm <- colnames(X)
  result <- susieR::univariate_regression(
    X, y,
    Z = Z, center = center, scale = scale,
    return_residuals = return_residuals
  )

  z_scores <- result$betahat / result$sebetahat
  p_values <- 2 * pnorm(abs(z_scores), lower.tail = FALSE)
  q_values <- compute_qvalues(p_values)

  result$betahat   <- setNames(result$betahat,   nm)
  result$sebetahat <- setNames(result$sebetahat, nm)
  result$z_scores  <- setNames(z_scores, nm)
  result$p_values  <- setNames(p_values, nm)
  result$q_values  <- setNames(q_values, nm)
  result
}


#' Perform Linear Regression for GWAS
#'
#' @param genotype Numeric matrix of genotypes (n x p), where n is the number of samples and p is the number of SNPs.
#' @param phenotype Numeric vector of phenotypes (length n).
#' @param covariates Optional numeric matrix of covariates (n x k), where k is the number of covariates.
#' @return A data frame containing regression results for each SNP, including \code{BETA}, \code{SE}, \code{Z}, \code{P}, and \code{Q}.
#' @examples
#' genotype <- matrix(rbinom(1000 * 10, 2, 0.3), 1000, 10)
#' phenotype <- rnorm(1000)
#' results <- run_linear_regression1(genotype, phenotype)
#' @noRd
run_linear_regression <- function(genotype, phenotype, covariates = NULL, phenotype_id = NULL) {
  if (!is.null(covariates)) {
    covariates <- as.data.frame(lapply(covariates, as.numeric))
  }

  reg_results <- univariate_regression(
    X = genotype,
    y = phenotype,
    Z = covariates,
    center = TRUE,
    scale = FALSE
  )

  snp_info <- parse_variant_id(colnames(genotype))

  data.frame(
    phenotype_id = if (!is.null(phenotype_id)) phenotype_id else NA,
    chr = snp_info$chrom,
    pos = snp_info$pos,
    A1 = snp_info$A1,
    A2 = snp_info$A2,
    variant_id = colnames(genotype),
    beta = reg_results$betahat,
    se = reg_results$sebetahat,
    z = reg_results$z_scores,
    p = reg_results$p_values,
    q = reg_results$q_values,
    N = nrow(genotype)
  )
}

#' Main QUAIL pipeline
#' QUAIL vQTL Analysis Pipeline
#'
#' @param genotype numeric matrix (n x p) of genotypes.
#' @param rank_score numeric vector (n x 1) of rank scores from Step 1.
#' @param phenotype optional numeric vector (n x 1) of original phenotype values.
#' @param covariates optional numeric matrix (n x k) of covariates.
#' @return A data frame containing vQTL results.
#' @export
#' @examples
#' \dontrun{
#' results <- QUAIL_pipeline(genotype, rank_score, covariates = covariates)
#' }
QUAIL_pipeline <- function(genotype, rank_score, phenotype = NULL,
                           covariates = NULL, phenotype_id = NULL) {
  start_time <- Sys.time()

  # Validate rank_score
  if (!is.numeric(rank_score)) {
    stop("rank_score must be a numeric vector.")
  }

  # Validate covariates
  if (!is.null(covariates)) {
    if (!is.data.frame(covariates)) {
      covariates <- as.data.frame(covariates)
    }
    covariates <- as.data.frame(lapply(covariates, as.numeric))
  }

  # Perform vQTL analysis
  vqtl_results <- run_linear_regression(genotype, rank_score, covariates, phenotype_id = phenotype_id)

  end_time <- Sys.time()
  cat("\nTotal vQTL runtime:", round(difftime(end_time, start_time, units = "secs"), 2), " seconds\n")

  return(vqtl_results)
}
