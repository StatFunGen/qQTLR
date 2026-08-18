#' @title Quantile TWAS Weight Calculation and QTL Analysis
#'
#' @description
#' This file contains functions for performing Quantile Transcriptome-Wide
#' Association Studies (TWAS) weight calculations and Quantile QTL analysis.
#' It provides tools for screening quantile regression results, performing
#' LD clumping and pruning, and calculating TWAS weights.
#'
#' @details
#' The main function in this file is `quantile_twas_weight_pipeline`, which
#' orchestrates the entire analysis process. Other functions are helper
#' functions used within the main pipeline.
#'

#' Qrank Score Test Screen
#' @param X Matrix of predictors
#' @param Y Matrix or vector of response variables
#' @param Z Matrix of covariates (optional)
#' @param tau.list Vector of quantiles to be analyzed
#' @param screen_threshold Significance threshold for adjusted p-values
#' @param screen_method Method for p-value adjustment ('fdr' or 'qvalue')
#' @param top_count Number of top SNPs to select
#' @param top_percent Percentage of top SNPs to select
#' @return A list containing various results from the QR screen
#' @importFrom tidyr separate
#' @importFrom dplyr %>% mutate select
#' @export
qr_screen <- function(
    X, Y, Z = NULL, tau.list = seq(0.05, 0.95, by = 0.05),
    screen_threshold = 0.05, screen_method = "qvalue", top_count = 10, top_percent = 15) {
  # Make sure quantreg is installed
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("To use this function, please install quantreg: https://cran.r-project.org/web/packages/quantreg/index.html")
  }
  p <- ncol(X)
  pvec <- rep(NA, p)
  ltau <- length(tau.list)
  quantile.pvalue <- matrix(NA, nrow = p, ncol = ltau, dimnames = list(colnames(X), paste("p_qr", tau.list, sep = "_")))
  quantile.zscore <- matrix(NA, nrow = p, ncol = ltau, dimnames = list(colnames(X), paste("zscore_qr", tau.list, sep = "_")))
  y <- as.matrix(Y)

  if (!is.null(Z)) {
    zz <- cbind(rep(1, nrow(y)), Z)
  } else {
    zz <- matrix(1, nrow = nrow(y), ncol = 1)
  }

  ranks_list <- lapply(tau.list, function(tau) {
    suppressWarnings(quantreg::rq.fit.br(zz, y, tau = tau)$dual - (1 - tau))
  })

  for (ip in 1:p) {
    x <- as.matrix(X[, ip])
    VN <- matrix(0, nrow = ltau, ncol = ltau)
    for (i in 1:ltau) {
      for (j in 1:ltau) {
        VN[i, j] <- min(tau.list[i], tau.list[j]) - tau.list[i] * tau.list[j]
      }
    }

    if (!is.null(Z)) {
      # xstar = lm(x ~ zz - 1)$residual
      xstar <- .lm.fit(zz, x)$residual
    } else {
      xstar <- x
    }

    SN <- NULL
    for (itau in 1:ltau) {
      Sn <- as.matrix(t(xstar) %*% ranks_list[[itau]])
      SN <- c(SN, Sn)
    }
    VN2 <- matrix(outer(VN, t(xstar) %*% xstar, "*"), nrow = ltau)
    z_score <- SN / sqrt(diag(VN2))
    pvalue1 <- pchisq(SN^2 / diag(VN2), 1, lower.tail = FALSE)
    names(pvalue1) <- tau.list
    quantile.pvalue[ip, ] <- pvalue1
    quantile.zscore[ip, ] <- z_score
    e <- solve(chol(VN2))
    SN2 <- t(e) %*% SN
    pvalue <- pchisq(sum(SN2^2), ltau, lower.tail = FALSE)
    pvec[ip] <- pvalue
  }

  pvec <- apply(quantile.pvalue, 1, pval_cauchy)

  if (screen_method == "fdr") {
    adjusted_pvalues <- p.adjust(pvec)
    method_col_name <- "fdr_p_qr"
    method_quantile_names <- paste0("fdr_p_qr_", tau.list)
    quantile_adjusted_pvalues <- apply(quantile.pvalue, 2, p.adjust)
  } else if (screen_method == "qvalue") {
    adjusted_pvalues <- compute_qvalues(pvec)
    method_col_name <- "qvalue_qr"
    method_quantile_names <- paste0("qvalue_qr_", tau.list)
    quantile_adjusted_pvalues <- apply(quantile.pvalue, 2, compute_qvalues)
  } else {
    stop("Invalid screen_method. Choose 'fdr' or 'qvalue'.")
  }

  # Ensure quantile_adjusted_pvalues is always a matrix (handles single-variant case)
  if (!is.matrix(quantile_adjusted_pvalues)) {
    quantile_adjusted_pvalues <- matrix(quantile_adjusted_pvalues, nrow = 1,
                                        dimnames = list(colnames(X), names(quantile_adjusted_pvalues)))
  }

  sig_SNP_threshold <- which(adjusted_pvalues < screen_threshold)
  sig_SNP_top_count <- order(adjusted_pvalues)[1:top_count]
  sig_SNP_top_percent <- order(adjusted_pvalues)[1:max(1, round(length(adjusted_pvalues) * top_percent / 100))]

  sig.SNPs_names <- colnames(X)[sig_SNP_threshold]
  sig.SNPs_names_top_count <- colnames(X)[sig_SNP_top_count]
  sig.SNPs_names_top_percent <- colnames(X)[sig_SNP_top_percent]
  phenotype_id <- colnames(y)[1]

  df_result <- data.frame(
    phenotype_id = phenotype_id,
    variant_id = colnames(X),
    p_qr = pvec
  )

  # Add quantile-specific p-values
  for (tau in tau.list) {
    df_result[[paste0("p_qr_", tau)]] <- quantile.pvalue[, paste0("p_qr_", tau)]
  }

  # Add overall q-value
  df_result[[method_col_name]] <- adjusted_pvalues

  # Add quantile-specific q-values
  for (tau in tau.list) {
    df_result[[paste0(method_col_name, "_", tau)]] <- quantile_adjusted_pvalues[, paste0("p_qr_", tau)]
  }

  # Add quantile-specific z-scores
  for (tau in tau.list) {
    df_result[[paste0("zscore_qr_", tau)]] <- quantile.zscore[, paste0("zscore_qr_", tau)]
  }

  # Split variant_id and reorder columns
  parsed <- parse_variant_id(df_result$variant_id)
  df_result <- df_result %>%
    mutate(chr = parsed$chrom, pos = parsed$pos, A2 = parsed$A2, A1 = parsed$A1)

  # Define the column order
  col_order <- c("chr", "pos", "A2", "A1", "phenotype_id", "variant_id", "p_qr")
  col_order <- c(col_order, paste0("p_qr_", tau.list))
  col_order <- c(col_order, method_col_name)
  col_order <- c(col_order, paste0(method_col_name, "_", tau.list))
  col_order <- c(col_order, paste0("zscore_qr_", tau.list))

  # Reorder the columns
  df_result <- df_result %>% select(all_of(col_order))

  return(list(
    df_result = df_result,
    tau_list = tau.list,
    quantile_pvalue = quantile.pvalue,
    quantile_zscore = quantile.zscore,
    pvec = pvec,
    adjusted_pvalues = adjusted_pvalues,
    sig_SNP_threshold = sig_SNP_threshold,
    sig.SNPs_names = sig.SNPs_names,
    sig_SNP_top_count = sig_SNP_top_count,
    sig_SNP_top_percent = sig_SNP_top_percent
  ))
}

#' Perform Clumping and Pruning Across Contexts (Quantiles)
#'
#' Two-stage LD clumping tailored to quantile-regression screens: first, a
#' per-tau p-value-weighted clump; second, a MAF-weighted (or uninformed)
#' final clump over the per-tau union. Delegates the underlying bigsnpr
#' call to \code{pecotmr::ld_clump_by_score}.
#'
#' @param X Matrix of genotypes (n x p_sig) restricted to significant SNPs.
#' @param qr_results Output of \code{qr_screen()}; must contain
#'   \code{sig.SNPs_names}, \code{tau_list}, \code{quantile_pvalue}, and
#'   \code{sig_SNP_threshold}.
#' @param maf_list Optional per-SNP minor allele frequencies aligned to the
#'   original (pre-screen) SNP ordering, used to score the final clump.
#' @param ld_clump_r2 r-squared threshold for per-tau clumping. Default 0.2.
#' @param final_clump_r2 r-squared threshold for the final MAF clump.
#'   Default 0.8.
#' @return A list with \code{final_SNPs} and \code{clumped_SNPs}.
#' @export
multicontext_ld_clumping <- function(X, qr_results, maf_list = NULL, ld_clump_r2 = 0.2, final_clump_r2 = 0.8) {
  sig_SNPs_names <- qr_results$sig.SNPs_names

  if (length(sig_SNPs_names) == 0) {
    return(list(final_SNPs = NULL, clumped_SNPs = NULL, message = "No significant SNPs found"))
  }

  if (ncol(X) == 1) {
    message("Only one SNP in X. Skipping LD clumping.")
    final_SNPs <- colnames(X)
    return(list(final_SNPs = final_SNPs, clumped_SNPs = final_SNPs,
                message = "Only one SNP, no LD clumping performed"))
  }

  # Parse chr/pos once.
  parsed_snp_info <- do.call(rbind, strsplit(sig_SNPs_names, ":"))
  chr <- as.numeric(sub("^chr", "", parsed_snp_info[, 1]))
  pos <- as.numeric(parsed_snp_info[, 2])

  # Stage 1: per-tau clump by -log10(p).
  clumped_snp_list <- list()
  for (itau in seq_along(qr_results$tau_list)) {
    tau_name <- paste0("p_qr_", qr_results$tau_list[itau])
    p_values_quantile <- qr_results$quantile_pvalue[, tau_name][qr_results$sig_SNP_threshold]
    log_p_values <- -log10(p_values_quantile)
    clumped_snp_list[[tau_name]] <- ld_clump_by_score(
      X, score = log_p_values, chr = chr, pos = pos, r2 = ld_clump_r2
    )
  }

  clumped_snp_union <- unique(unlist(clumped_snp_list))
  clumped_SNPs_name <- sig_SNPs_names[clumped_snp_union]
  message("Number of SNPs after union of clumping: ", length(clumped_snp_union))

  if (length(clumped_snp_union) == 1) {
    message("Only one SNP found in the union. Skipping final LD clumping.")
    final_SNPs <- sig_SNPs_names[clumped_snp_union]
    return(list(final_SNPs = final_SNPs, clumped_SNPs = clumped_SNPs_name))
  }

  # Stage 2: final clump over the union, ordered by chr/pos, scored by MAF (or NULL).
  sorted_indices <- order(chr[clumped_snp_union], pos[clumped_snp_union])
  chr_sorted <- chr[clumped_snp_union][sorted_indices]
  pos_sorted <- pos[clumped_snp_union][sorted_indices]

  maf_values <- NULL
  if (!is.null(maf_list)) {
    maf_values <- maf_list[qr_results$sig_SNP_threshold][clumped_snp_union][sorted_indices]
  }
  G_union <- X[, clumped_snp_union, drop = FALSE][, sorted_indices]

  final_clumped <- ld_clump_by_score(
    G_union, score = maf_values, chr = chr_sorted, pos = pos_sorted, r2 = final_clump_r2
  )

  final_SNPs <- sig_SNPs_names[clumped_snp_union][sorted_indices][final_clumped]
  message("Number of final SNPs after MAF-based clumping (if applied): ", length(final_SNPs))

  list(final_SNPs = final_SNPs, clumped_SNPs = clumped_SNPs_name)
}

#' Perform Quantile Regression Analysis to get beta
#' @param X Matrix of predictors
#' @param Y Matrix or vector of response variables
#' @param Z Matrix of covariates (optional)
#' @param tau_values Vector of quantiles to be analyzed
#' @return A data frame with QR coefficients for each quantile
#' @importFrom tidyr pivot_wider separate
#' @importFrom dplyr %>% mutate select
#' @export
perform_qr_analysis <- function(X, Y, Z = NULL, tau_values = seq(0.05, 0.95, by = 0.05)) {
  # Make sure quantreg is installed
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("To use this function, please install quantreg: https://cran.r-project.org/web/packages/quantreg/index.html")
  }
  # Convert Y and X to matrices if they aren't already
  pheno.mat <- as.matrix(Y)
  geno.mat <- as.matrix(X)

  # Initialize an empty result table to store results
  result_table <- data.frame(
    phenotype_id = character(),
    variant_id = character(),
    tau = numeric(),
    predictor_coef = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over each tau value to perform quantile regression
  for (tau in tau_values) {
    # Loop over each SNP/variant in geno.mat (X)
    for (n in 1:ncol(geno.mat)) {
      response <- pheno.mat # Y
      predictor <- geno.mat[, n] # X
      phenotype_id <- colnames(pheno.mat)
      variant_id <- colnames(geno.mat)[n]

      # Construct the design matrix based on whether Z is provided
      if (is.null(Z)) {
        # If no covariates, include intercept and predictor
        X_design <- cbind(1, predictor)
      } else {
        # If covariates are provided, include them in the design matrix
        X_design <- cbind(1, predictor, as.matrix(Z))
      }

      # Fit the quantile regression model using rq.fit.br
      mod <- suppressWarnings(quantreg::rq.fit.br(X_design, response, tau = tau))

      # Extract the coefficient for the predictor (second coefficient)
      predictor_coef <- mod$coefficients[2] # Coefficient for predictor

      # Create a row with the results and append to the result table
      row <- data.frame(
        phenotype_id = phenotype_id,
        variant_id = variant_id,
        tau = tau,
        predictor_coef = predictor_coef,
        stringsAsFactors = FALSE
      )
      result_table <- rbind(result_table, row)
    }
  }

  # Reshape result_table to a wide format, so each tau's results are in separate columns
  result_table_wide <- result_table %>%
    pivot_wider(
      id_cols = c(phenotype_id, variant_id),
      names_from = tau,
      values_from = predictor_coef,
      names_prefix = "coef_qr_"
    )
  parsed_ids <- parse_variant_id(result_table_wide$variant_id)
  result_table_wide <- result_table_wide %>%
    mutate(chr = parsed_ids$chrom, pos = parsed_ids$pos, A2 = parsed_ids$A2, A1 = parsed_ids$A1) %>%
    select(chr, pos, A2, A1, everything())

  # Return the wide format result table
  return(result_table_wide)
}


#' Calculate QR Coefficients and Pseudo R-squared Across Multiple Quantiles
#'
#' This function calculates quantile regression coefficients and pseudo R-squared values across multiple quantiles,
#' while handling problematic columns that might affect the rank of the design matrix.
#'
#' @param AssocData List containing X, Y, C, and X.filter
#' @param tau.list Vector of quantiles to be analyzed
#' @param strategy The strategy for removing problematic columns ("variance", "correlation", or "response_correlation")
#' @return A list containing the cleaned X matrix, beta matrix as twas weight, and pseudo R-squared values
#' @noRd
calculate_qr_and_pseudo_R2 <- function(AssocData, tau.list, strategy = c("correlation", "variance", "response_correlation"),
                                       corr_thresholds = seq(0.75, 0.5, by = -0.05)) {
  # Make sure quantreg is installed
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("To use this function, please install quantreg: https://cran.r-project.org/web/packages/quantreg/index.html")
  }
  strategy <- match.arg(strategy)
  # Check and handle problematic columns affecting the full rank of the design matrix
  AssocData$X.filter <- enforce_design_full_rank(X = AssocData$X.filter, C = AssocData$C, strategy = strategy, response = AssocData$Y, corr_thresholds = corr_thresholds)
  snp_names <- colnames(AssocData$X.filter)
  # Build the cleaned design matrix using the filtered X and unnamed C

  # Fit the models for all tau values
  message("Start fitting full model for all taus...")
  fit_full <- suppressWarnings(quantreg::rq(Y ~ X.filter + C, tau = tau.list, data = AssocData))
  message("Finished fitting full model. Start fitting intercept-only model for all taus...")
  fit_intercept <- suppressWarnings(quantreg::rq(AssocData$Y ~ 1, tau = tau.list, data = AssocData))
  message("Finished fitting intercept-only model.")
  # Define the rho function for pseudo R2 calculation
  rho <- function(u, tau) {
    u * (tau - (u < 0))
  }

  # Prepare to store the pseudo R2 results
  pseudo_R2 <- numeric(length(tau.list))
  names(pseudo_R2) <- tau.list

  # Calculate pseudo R2 for each tau
  for (i in seq_along(tau.list)) {
    tau <- tau.list[i]

    # Get residuals for the intercept-only and full models
    residuals0 <- residuals(fit_intercept, subset = i)
    residuals1 <- residuals(fit_full, subset = i)

    # Calculate and store pseudo R2 for each tau
    rho0 <- sum(rho(residuals0, tau))
    rho1 <- sum(rho(residuals1, tau))
    pseudo_R2[i] <- 1 - rho1 / rho0
  }

  # Extract the coefficients for the SNPs
  num_filter_vars <- ncol(AssocData$X.filter)
  beta_mat <- coef(fit_full)[2:(1 + num_filter_vars), , drop = FALSE]
  rownames_beta <- rownames(beta_mat)
  if (ncol(AssocData$X.filter) == 1) {
    rownames(beta_mat) <- snp_names
  } else {
    rownames_beta <- rownames(beta_mat)
    rownames(beta_mat) <- gsub("^X.filter", "", rownames_beta)
  }
  return(list(X.filter = AssocData$X.filter, beta_mat = beta_mat, pseudo_R2 = pseudo_R2))
}

#' Calculate Heterogeneity of Beta Coefficients Across Quantiles
#'
#' This function calculates the heterogeneity of beta coefficients across multiple quantiles for each variant_id.
#' Heterogeneity is computed as log(sd(beta) / abs(mean(beta))).
#'
#' @param rq_coef_result Data frame containing variant_id and QR coefficient columns
#' @return A data frame with variant_id and heterogeneity values
#' @noRd
calculate_coef_heterogeneity <- function(rq_coef_result) {
  # Identify all the columns starting with "coef_qr_" (quantile regression coefficient columns)
  coef_cols <- grep("^coef_qr_", colnames(rq_coef_result), value = TRUE)

  # Create a new data frame with variant_id and heterogeneity
  heterogeneity_result <- data.frame(
    variant_id = rq_coef_result$variant_id,
    coef_heter = apply(rq_coef_result[, coef_cols], 1, function(beta) {
      # Compute the mean and standard deviation, ignoring NAs
      beta_mean <- mean(beta, na.rm = TRUE)
      beta_sd <- sd(beta, na.rm = TRUE)

      # Handle the case where mean(beta) is 0 to avoid division by zero
      if (abs(beta_mean) == 0) {
        return(NA) # Return NA if mean is zero
      }

      # Compute the heterogeneity: log(sd(beta) / abs(mean(beta)))
      heterogeneity <- log(beta_sd / abs(beta_mean))
      return(heterogeneity)
    }),
    stringsAsFactors = FALSE
  )

  # Return only variant_id and heterogeneity
  return(heterogeneity_result)
}

#' Calculate Xi Correlation for QR Coefficients
#'
#' This function calculates the xi correlation coefficient and p-value for each variant,
#' measuring the functional dependence between tau values and QR coefficients.
#' Uses coef_qr_0.1 to coef_qr_0.9 (17 values, excluding 0.05 and 0.95).
#'
#' @param rq_coef_result Data frame containing variant_id and coef_qr_* columns
#' @param tau_range Numeric vector of tau values to use (default: seq(0.1, 0.9, by = 0.05))
#' @param min_valid Minimum number of valid (non-NA) coefficients required (default: 10)
#' @return A data frame with variant_id, xi, and xi_pval columns
#' @export
calculate_xi_correlation <- function(rq_coef_result, tau_range = seq(0.1, 0.9, by = 0.05), min_valid = 10) {
  if (!requireNamespace("XICOR", quietly = TRUE)) {
    stop("Package 'XICOR' is required for xi correlation calculation. Please install it.")
  }
  # Build column names for the specified tau range
  coef_col_names <- paste0("coef_qr_", tau_range)

  # Check which columns exist
  existing_cols <- coef_col_names[coef_col_names %in% colnames(rq_coef_result)]

  if (length(existing_cols) == 0) {
    warning("No coef_qr columns found in the specified tau range")
    return(data.frame(
      variant_id = rq_coef_result$variant_id,
      xi = NA_real_,
      xi_pval = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  # Extract tau values from existing columns
  existing_tau <- as.numeric(gsub("coef_qr_", "", existing_cols))

  # Calculate xi for each variant
  xi_results <- apply(rq_coef_result[, existing_cols, drop = FALSE], 1, function(coef_values) {
    # Get valid (non-NA) coefficients
    valid_indices <- !is.na(coef_values)
    valid_coefs <- as.numeric(coef_values[valid_indices])
    valid_tau <- existing_tau[valid_indices]

    # Check if enough valid values
    if (length(valid_coefs) < min_valid) {
      return(c(xi = NA_real_, xi_pval = NA_real_))
    }

    # Calculate xi correlation
    tryCatch({
      xicor_result <- XICOR::xicor(valid_tau, y = valid_coefs, pvalue = TRUE, method = "asymptotic")
      return(c(xi = xicor_result$xi, xi_pval = xicor_result$pval))
    }, error = function(e) {
      return(c(xi = NA_real_, xi_pval = NA_real_))
    })
  })

  # Convert to data frame
  xi_df <- data.frame(
    variant_id = rq_coef_result$variant_id,
    xi = xi_results["xi", ],
    xi_pval = xi_results["xi_pval", ],
    stringsAsFactors = FALSE
  )

  return(xi_df)
}

#' Bootstrap-Calibrated p-value for the Xi Correlation Heterogeneity Test
#'
#' \code{calculate_xi_correlation}'s asymptotic p-value
#' (\code{XICOR::xicor(..., method = "asymptotic")}) is badly miscalibrated for this use
#' case: the tau-grid coefficients \eqn{\hat\beta(\tau)} are all estimated from the SAME n
#' individuals, so adjacent-tau estimates are highly correlated. Chatterjee's xi is built
#' for i.i.d. (X,Y) pairs and mistakes this sampling-induced smoothness for genuine
#' tau-dependence: in simulation under a true null (constant per-allele effect, no real
#' heterogeneity), the asymptotic p-value rejects at nominal alpha = 0.05 in roughly 87-90\%
#' of replicates (see \code{manuscript/heter_test_khmaladze/task1_calibrated_xi_bootstrap.R}
#' for the verification). Naively permuting the tau <-> beta_hat(tau) correspondence does
#' NOT fix this (still ~87\% false-positive rate) -- permuting labels on one already
#' correlated, fixed vector never reconstructs a valid null sampling distribution.
#'
#' This function instead regenerates the null and refits the FULL quantile regression
#' process on each bootstrap replicate: (1) fit a null-constrained (median location-shift)
#' model to get a constant-effect fit and residuals, (2) resample residuals with replacement
#' and add them back onto the null fitted values (same fixed genotype/covariate design),
#' (3) refit all tau levels on the synthetic outcome and recompute xi, (4) repeat n_boot
#' times to build an empirical null distribution specific to this variant's actual sample.
#'
#' Because of this refit-based design, calibration CANNOT be computed post-hoc from
#' already-saved \code{coef_qr_*} columns alone (unlike \code{calculate_xi_correlation}) --
#' it needs the individual-level \code{X}, \code{Y}, \code{Z}. It is also \code{O(n_boot)}
#' times more expensive per variant than the asymptotic version, so it should only be
#' requested for a restricted candidate set (e.g. via \code{keep_variants} upstream), not
#' genome-wide.
#'
#' @param X Genotype matrix, one column per variant, column names = variant_id
#' @param Y Single-column phenotype matrix/vector
#' @param Z Covariate matrix (optional)
#' @param taus Quantile grid (default: the pipeline's 19-point grid)
#' @param n_boot Number of null-resampling bootstrap replicates (default 199)
#' @return A data frame with variant_id, xi_obs, xi_pval_asymp (for comparison against the
#'   uncalibrated version), and xi_pval_calibrated
#' @export
calculate_xi_correlation_calibrated <- function(X, Y, Z = NULL,
                                                  taus = seq(0.05, 0.95, by = 0.05),
                                                  n_boot = 199) {
  if (!requireNamespace("XICOR", quietly = TRUE)) {
    stop("Package 'XICOR' is required for xi correlation calculation. Please install it.")
  }
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for xi correlation calculation. Please install it.")
  }
  geno.mat <- as.matrix(X)
  y_vec <- as.numeric(as.matrix(Y)[, 1])
  variant_ids <- colnames(geno.mat)
  z_mat <- if (is.null(Z)) NULL else as.matrix(Z)

  design_for <- function(g_col) if (is.null(z_mat)) cbind(1, g_col) else cbind(1, g_col, z_mat)

  fit_betas_g <- function(g_col, y, taus) {
    Xd <- design_for(g_col)
    sapply(taus, function(tt) {
      mod <- suppressWarnings(quantreg::rq.fit.br(Xd, y, tau = tt))
      mod$coefficients[2]
    })
  }

  message("Computing bootstrap-calibrated xi p-values (n_boot=", n_boot,
          "); this refits the full ", length(taus),
          "-tau quantile process n_boot times PER VARIANT -- restrict to a small candidate set.")

  res <- lapply(seq_len(ncol(geno.mat)), function(j) {
    g_col <- geno.mat[, j]
    Xd <- design_for(g_col)
    med_mod <- suppressWarnings(quantreg::rq.fit.br(Xd, y_vec, tau = 0.5))
    fitted_null <- as.numeric(Xd %*% med_mod$coefficients)
    resid0 <- y_vec - fitted_null

    betas_obs <- fit_betas_g(g_col, y_vec, taus)
    obs_fit <- XICOR::xicor(taus, betas_obs, method = "asymptotic", pvalue = TRUE)
    xi_obs <- obs_fit$xi

    xi_boot <- sapply(seq_len(n_boot), function(b) {
      y_star <- fitted_null + sample(resid0, length(resid0), replace = TRUE)
      betas_star <- fit_betas_g(g_col, y_star, taus)
      XICOR::xicor(taus, betas_star, method = "asymptotic")
    })
    c(xi_obs = unname(xi_obs), xi_pval_asymp = unname(obs_fit$pval),
      xi_pval_calibrated = unname((1 + sum(xi_boot >= xi_obs)) / (1 + n_boot)))
  })

  res_mat <- do.call(rbind, res)
  data.frame(variant_id = variant_ids, res_mat, stringsAsFactors = FALSE)
}

#' KhmaladzeTest Heterogeneity Statistic for Quantile Regression Coefficients
#'
#' Tests whether the genotype's quantile regression slope is constant across tau, using
#' the Khmaladze (1981) martingale-transform test as adapted to quantile regression by
#' Koenker \& Xiao (2002), via \code{quantreg::KhmaladzeTest()}. In simulation this test is
#' much better calibrated under a true null than Chatterjee's xi (~5\% rejection at nominal
#' alpha=0.05 vs. ~87-90\% for xi's asymptotic p-value), with good power against both
#' monotonic and non-monotonic ("hump") true heterogeneity, and reduced (but non-zero) power
#' for weak, tail-localized signal (see
#' \code{xi_calibration_check/xi_calibration_verification.md}, section 8, for the
#' verification this function's design is based on).
#'
#' Unlike \code{calculate_xi_correlation}, this CANNOT be computed from already-fitted
#' \code{coef_qr_*} columns -- \code{KhmaladzeTest} refits the full quantile regression
#' PROCESS from individual-level data (\code{rqProcess()} plus an \code{akj()} nonparametric
#' conditional-density estimate internally), so \code{X}, \code{Y}, \code{Z} must be supplied.
#'
#' Handling of covariates \code{Z} (revised 2026-08-10, credit Roulan Jiang for the
#' diagnosis): \code{Z} is now put DIRECTLY into the \code{KhmaladzeTest} formula
#' (\code{Y ~ G + Z1 + ... + Zp}), matching \code{calculate_xi_correlation_calibrated}'s
#' design instead of diverging from it. An earlier version of this function partialled
#' \code{Z} out of \code{Y} and \code{G} via OLS first (Frisch-Waugh-Lovell) because the
#' direct formula threw "system is computationally singular" on real data (~40-90
#' covariates) -- but that turned out to be caused by \code{summary.rq}'s default sparsity
#' estimator (\code{se="nid"}, which \code{KhmaladzeTest}'s \code{...} passes through),
#' not by the covariate count itself. Passing \code{se="ker"} instead avoids the crash
#' entirely (confirmed on TMEM106B/chr7:12086129:A:G, ROSMAP Exc_mega, n=737, 65
#' covariates: \code{se="nid"} still throws the singular error, \code{se="ker"} returns
#' \code{Tn=58.3} cleanly) and also reduces MAF-driven Type-I-error inflation under a true
#' null (MAF=0.05: ~10\% down to ~6\% at nominal alpha=0.05 in simulation, vs. \code{"nid"}'s
#' ~15\%). \code{se="ker"} is now the default here (can be overridden via \code{...}).
#'
#' Because \code{Z} is back in the formula, \code{Tn} (joint, all slopes) and \code{THn}
#' (named per-slope vector) are no longer numerically identical -- use \code{THn["G"]}
#' (returned as \code{khm_THn_G}), NOT \code{Tn}, exactly as originally documented: \code{Tn}
#' would reject whenever ANY covariate has its own tau-varying effect, regardless of whether
#' the genotype term itself is heterogeneous.
#'
#' No exact analytic p-value for \code{THn} was locatable in this environment (the tabulated
#' critical values from Koenker \& Xiao 2002 were not available); set \code{n_boot > 0} to
#' obtain an empirical, null-consistent p-value via the same residual-bootstrap-and-refit
#' principle as \code{calculate_xi_correlation_calibrated} -- the null-generation logic is now
#' IDENTICAL between the two functions (fit a median/location-shift model on the FULL design
#' including \code{Z}, resample residuals onto the same fixed G/Z, refit all tau, recompute
#' the statistic), differing only in which statistic (xi vs. THn["G"]) is recomputed each
#' replicate. This is expensive (\code{n_boot} full \code{KhmaladzeTest} refits per variant,
#' each itself more costly than a plain \code{rq} fit) -- only enable it for a small
#' candidate set, not genome-wide.
#'
#' KNOWN OPEN LIMITATION (raised by Roulan Jiang, unresolved as of 2026-08-10): the
#' null-generating model for the bootstrap treats \code{Z}'s coefficients as FIXED at
#' their median-regression values, i.e. it assumes \code{Z} itself has no tau-varying
#' effect either. If a covariate genuinely does have a tau-varying effect, the synthetic
#' null data does not reflect that, and it is not yet established whether this biases
#' \code{khm_pval_calibrated} (or \code{xi_pval_calibrated}) upward or downward. Roulan is
#' investigating an alternative (Chernozhukov \& Fernandez-Val 2005 subsampling inference
#' for the QR process) that may not share this assumption; not yet implemented here.
#'
#' @param X Genotype matrix, one column per variant, column names = variant_id
#' @param Y Single-column phenotype matrix/vector
#' @param Z Covariate matrix (optional)
#' @param taus Quantile grid (default: the pipeline's 19-point grid)
#' @param n_boot Number of null-resampling bootstrap replicates for an empirical p-value per
#'   variant (default 0 = skip, return only the Tn/THn_G statistics)
#' @param se Sparsity/density estimator passed to \code{KhmaladzeTest} (default \code{"ker"};
#'   quantreg's default \code{"nid"} is prone to "computationally singular" errors with many
#'   covariates and somewhat worse Type-I-error control at low MAF -- see Details)
#' @return A data frame with variant_id, khm_Tn (joint statistic, for reference only),
#'   khm_THn_G (genotype-specific statistic -- use this one), and khm_pval_calibrated
#'   (only if n_boot > 0)
#' @export
calculate_khmaladze_heterogeneity <- function(X, Y, Z = NULL,
                                               taus = seq(0.05, 0.95, by = 0.05),
                                               n_boot = 0, se = "ker") {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for the KhmaladzeTest heterogeneity test. Please install it.")
  }
  geno.mat <- as.matrix(X)
  y_vec <- as.numeric(as.matrix(Y)[, 1])
  variant_ids <- colnames(geno.mat)
  z_mat <- if (is.null(Z)) NULL else as.matrix(Z)
  if (!is.null(z_mat)) colnames(z_mat) <- paste0("Zcov", seq_len(ncol(z_mat)))

  khm_formula <- if (is.null(z_mat)) {
    Y ~ G
  } else {
    stats::as.formula(paste("Y ~ G +", paste(colnames(z_mat), collapse = " + ")))
  }

  fit_one <- function(y, g_col) {
    df <- data.frame(Y = y, G = g_col)
    if (!is.null(z_mat)) df <- cbind(df, as.data.frame(z_mat))
    T <- tryCatch(quantreg::KhmaladzeTest(khm_formula, data = df, taus = taus, nullH = "location", se = se),
                  error = function(e) NULL)
    if (is.null(T)) return(c(Tn = NA_real_, THn_G = NA_real_))
    thn_g <- if (is.null(names(T$THn))) T$Tn else unname(T$THn["G"])  # single-predictor case: THn is unnamed
    c(Tn = unname(T$Tn), THn_G = thn_g)
  }

  khm_stats <- t(sapply(seq_len(ncol(geno.mat)), function(j) fit_one(y_vec, geno.mat[, j])))
  out <- data.frame(variant_id = variant_ids,
                     khm_Tn = khm_stats[, "Tn"],
                     khm_THn_G = khm_stats[, "THn_G"],
                     stringsAsFactors = FALSE)

  if (n_boot > 0) {
    message("Computing bootstrap-calibrated KhmaladzeTest p-values (n_boot=", n_boot,
            "); this refits the full quantile process n_boot times PER VARIANT -- restrict to a small candidate set.")
    design_for <- function(g_col) if (is.null(z_mat)) cbind(1, g_col) else cbind(1, g_col, z_mat)
    out$khm_pval_calibrated <- sapply(seq_len(ncol(geno.mat)), function(j) {
      obs_THn <- out$khm_THn_G[j]
      if (is.na(obs_THn)) return(NA_real_)
      g_col <- geno.mat[, j]
      Xd <- design_for(g_col)
      # null-constrained (median location-shift) model on the FULL design (G and Z), same
      # null-generation logic as calculate_xi_correlation_calibrated
      med_mod <- suppressWarnings(quantreg::rq.fit.br(Xd, y_vec, tau = 0.5))
      fitted_null <- as.numeric(Xd %*% med_mod$coefficients)
      resid0 <- y_vec - fitted_null
      boot_THn <- sapply(seq_len(n_boot), function(b) {
        y_star <- fitted_null + sample(resid0, length(resid0), replace = TRUE)
        fit_one(y_star, g_col)["THn_G"]
      })
      (1 + sum(boot_THn >= obs_THn, na.rm = TRUE)) / (1 + sum(!is.na(boot_THn)))
    })
  }
  out
}

#' Wald Test for Equality of Quantile Regression Slopes Across Tau
#'
#' Tests H0: beta(tau_1) = ... = beta(tau_L) via \code{quantreg::anova.rq(fit, joint =
#' FALSE)} (Bassett \& Koenker 1982; Koenker \& Bassett 1982) -- suggested by Tianying Wang.
#' Unlike \code{calculate_xi_correlation}/\code{calculate_xi_correlation_calibrated}, this
#' statistic's asymptotic covariance estimator is built for the STACKED beta_hat(tau) vector
#' from the start, so it does not ignore the across-tau correlation the way xi's asymptotic
#' formula does -- \strong{no bootstrap is needed}. In simulation (n=700, MAF=0.3, same 4
#' scenarios as \code{calculate_xi_correlation_calibrated}'s verification): true-null
#' rejection ~3-4\% at alpha=0.05 (well calibrated), ~100\% power for both monotonic and
#' non-monotonic ("hump") true heterogeneity, ~13\% power for a weak one-sided (qQTL-only-
#' shaped) signal -- calibration and power at least as good as the bootstrap methods, at a
#' small fraction of the computational cost (a single \code{rq} fit + \code{anova} call per
#' variant, no refitting).
#'
#' MAF sensitivity (open caveat, mirrors \code{calculate_khmaladze_heterogeneity}): under a
#' true null, rejection rate at alpha=0.05 rises as MAF falls (MAF=0.05: ~6-10\%; MAF=0.5:
#' 0-2\%), because \code{anova.rq}'s sparsity/density estimator is less stable for the
#' minor-allele group when it is small. Passing \code{se="ker"} (the default here) reduces
#' this inflation substantially versus quantreg's own default \code{se="nid"} (confirmed in
#' simulation: MAF=0.05 true-null rejection ~10\% with "nid" vs. ~6\% with "ker"; "nid" also
#' throws "computationally singular" on real data with many covariates in some cases, which
#' "ker" avoids -- credit Roulan Jiang for this diagnosis). Even with "ker", low-MAF variants
#' should be interpreted with some caution; this is not fully resolved.
#'
#' The p-value for the genotype column is read off \code{anova_res$table} by ROW NAME
#' (\code{"G"}), not by position -- indexing by position (e.g. \code{table$pvalue[1]}) is
#' only safe if the genotype term is guaranteed to be the first slope in the formula.
#'
#' @param X Genotype matrix, one column per variant, column names = variant_id
#' @param Y Single-column phenotype matrix/vector
#' @param Z Covariate matrix (optional) -- included directly in the \code{rq} formula (no
#'   Frisch-Waugh-Lovell residualization needed; confirmed not to crash on real data with
#'   ~90 covariates, n=737, unlike the earlier issue found with \code{KhmaladzeTest})
#' @param taus Quantile grid (default: the pipeline's 19-point grid)
#' @param se Sparsity/density estimator passed to \code{anova.rq} (default \code{"ker"})
#' @return A data frame with variant_id, wald_Tn (test statistic), wald_ndf/wald_ddf
#'   (numerator/denominator degrees of freedom), and wald_pval
#' @export
calculate_wald_heterogeneity <- function(X, Y, Z = NULL,
                                          taus = seq(0.05, 0.95, by = 0.05),
                                          se = "ker") {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for the Wald heterogeneity test. Please install it.")
  }
  geno.mat <- as.matrix(X)
  y_vec <- as.numeric(as.matrix(Y)[, 1])
  variant_ids <- colnames(geno.mat)
  z_mat <- if (is.null(Z)) NULL else as.matrix(Z)
  if (!is.null(z_mat)) colnames(z_mat) <- paste0("Zcov", seq_len(ncol(z_mat)))

  wald_formula <- if (is.null(z_mat)) {
    Y ~ G
  } else {
    stats::as.formula(paste("Y ~ G +", paste(colnames(z_mat), collapse = " + ")))
  }

  res <- t(sapply(seq_len(ncol(geno.mat)), function(j) {
    df <- data.frame(Y = y_vec, G = geno.mat[, j])
    if (!is.null(z_mat)) df <- cbind(df, as.data.frame(z_mat))
    fit <- tryCatch(suppressWarnings(quantreg::rq(wald_formula, data = df, tau = taus)),
                     error = function(e) NULL)
    a <- if (is.null(fit)) NULL else tryCatch(suppressWarnings(stats::anova(fit, joint = FALSE, se = se)),
                                                error = function(e) NULL)
    if (is.null(a)) return(c(Tn = NA_real_, ndf = NA_real_, ddf = NA_real_, pval = NA_real_))
    c(Tn = a$table["G", "Tn"], ndf = a$table["G", "ndf"], ddf = a$table["G", "ddf"], pval = a$table["G", "pvalue"])
  }))

  data.frame(variant_id = variant_ids,
             wald_Tn = res[, "Tn"], wald_ndf = res[, "ndf"], wald_ddf = res[, "ddf"],
             wald_pval = res[, "pval"], stringsAsFactors = FALSE)
}

#' Subsampling Test for Equality of Quantile Regression Slopes Across Tau
#'
#' Tests H0: beta(tau_1) = ... = beta(tau_L) using the recentered ("mimicking") subsampling
#' procedure of Chernozhukov \& Fernandez-Val (2005, \emph{Subsampling Inference on Quantile
#' Regression Processes}, Sankhya 67(2):253-276, their "Example 2: constant vs. heterogeneous
#' effects"), suggested by Roulan Jiang as an alternative to \code{calculate_khmaladze_heterogeneity}/
#' \code{calculate_wald_heterogeneity} that avoids density/sparsity estimation altogether (no
#' \code{se="nid"}/\code{se="ker"} choice, hence no analog of the crashes or MAF-driven Type-I-error
#' inflation those two functions document).
#'
#' \strong{Procedure}: let \eqn{\hat\beta_G(\tau)} be the genotype slope from a single
#' quantile regression process fit on the FULL sample (\code{Y ~ G + Z1 + ... + Zp}, all
#' taus). The observed inference process is \eqn{v_n(\tau) = \hat\beta_G(\tau) - \bar r}, where
#' \eqn{\bar r} is the grid average of \eqn{\hat\beta_G(\tau)} (the paper's estimate of the
#' common effect under H0). Then, \code{n_boot} times: draw a subsample of size
#' \code{b = round(subsample_frac * n)} \emph{without replacement} from the n individuals,
#' refit the SAME quantile process on just those b individuals, and take
#' \eqn{d_i(\tau) = \hat\beta_G^{(i)}(\tau) - \hat\beta_G(\tau)} -- the deviation of the
#' subsample's curve from the FULL-sample curve (this differencing against the full-sample
#' curve, rather than against zero or against the subsample's own null, is the paper's
#' "recentering"/"mimicking" trick: it makes the resampled statistic's distribution
#' consistently estimate the true null distribution of the test statistic under both the null
#' and local alternatives, which is what gives this test better finite-sample power than
#' textbook/canonical subsampling and makes it insensitive to the exact choice of
#' \code{subsample_frac}). Pointwise variances \eqn{Var_\tau} of \eqn{d_i(\tau)} across the
#' \code{n_boot} replicates are used as an Anderson-Darling-type weight (the paper's own
#' preferred choice, Section 3.3) -- this is what replaces the sparsity/density estimation
#' that \code{KhmaladzeTest}/\code{anova.rq} need. The test statistic (\code{statistic =
#' "cms"}, Cramer-von-Mises-Smirnov type, default) is
#' \eqn{S_n = (n/b) \sum_\tau v_n(\tau)^2 / Var_\tau(\tau)} (grid-averaged), or
#' (\code{statistic = "ks"}, Kolmogorov-Smirnov type) \eqn{S_n = \sqrt{n/b} \max_\tau
#' |v_n(\tau)| / \sqrt{Var_\tau(\tau)}}; the subsampling replicates give the SAME statistic
#' computed on \eqn{d_i(\tau)} in place of \eqn{v_n(\tau)}, and the p-value is the rank of
#' \eqn{S_n} among the \code{n_boot} replicate values.
#'
#' \strong{Relation to the open covariate critique (Roulan Jiang, raised for
#' \code{calculate_xi_correlation_calibrated}/\code{calculate_khmaladze_heterogeneity}):} those
#' two functions generate their null by resampling residuals onto a model that fixes \code{Z}'s
#' coefficients at their median-regression (tau=0.5) values, which implicitly assumes \code{Z}
#' itself has no real tau-varying effect. This function makes no such assumption -- subsamples
#' are drawn directly from the observed data with no synthetic null construction, so any real
#' tau-varying effect that \code{Z} has is preserved exactly the same way in both the
#' full-sample statistic and every subsample replicate. This is the main reason Roulan
#' suggested investigating it as a candidate replacement for the residual-bootstrap approach.
#'
#' Not yet cross-validated as extensively as the other three options here -- see
#' \code{manuscript/heter_test_khmaladze/task27_subsampling_verification.R} for the
#' Type-I-error/power simulation this implementation was checked against before use on real
#' data.
#'
#' @param X Genotype matrix, one column per variant, column names = variant_id
#' @param Y Single-column phenotype matrix/vector
#' @param Z Covariate matrix (optional)
#' @param taus Quantile grid (default: the pipeline's 19-point grid)
#' @param n_boot Number of subsample replicates (default 199)
#' @param subsample_frac Subsample size as a fraction of n (default 0.5); per the paper's
#'   recentering result this choice should not matter much for the test's size/power
#' @param statistic \code{"cms"} (Cramer-von-Mises-Smirnov, default, aggregates the whole
#'   curve) or \code{"ks"} (Kolmogorov-Smirnov, driven by the single most extreme tau)
#' @return A data frame with variant_id, subsamp_Sn (observed statistic), subsamp_b
#'   (subsample size actually used), and subsamp_pval
#' @export
calculate_subsampling_heterogeneity <- function(X, Y, Z = NULL,
                                                  taus = seq(0.05, 0.95, by = 0.05),
                                                  n_boot = 199, subsample_frac = 0.5,
                                                  statistic = c("cms", "ks")) {
  statistic <- match.arg(statistic)
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("Package 'quantreg' is required for the subsampling heterogeneity test. Please install it.")
  }
  geno.mat <- as.matrix(X)
  y_vec <- as.numeric(as.matrix(Y)[, 1])
  n <- length(y_vec)
  variant_ids <- colnames(geno.mat)
  z_mat <- if (is.null(Z)) NULL else as.matrix(Z)
  if (!is.null(z_mat)) colnames(z_mat) <- paste0("Zcov", seq_len(ncol(z_mat)))
  b <- max(round(subsample_frac * n), (if (is.null(z_mat)) 0 else ncol(z_mat)) + 5)

  ss_formula <- if (is.null(z_mat)) {
    Y ~ G
  } else {
    stats::as.formula(paste("Y ~ G +", paste(colnames(z_mat), collapse = " + ")))
  }

  beta_g_curve <- function(idx, g_col) {
    df <- data.frame(Y = y_vec[idx], G = g_col[idx])
    if (!is.null(z_mat)) df <- cbind(df, as.data.frame(z_mat[idx, , drop = FALSE]))
    fit <- tryCatch(suppressWarnings(quantreg::rq(ss_formula, data = df, tau = taus)),
                     error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, length(taus)))
    cf <- stats::coef(fit)
    if (is.matrix(cf)) unname(cf["G", ]) else unname(cf["G"])
  }

  message("Computing subsampling-calibrated heterogeneity p-values (n_boot=", n_boot,
          ", subsample size b=", b, " of n=", n,
          "); this refits the full quantile process n_boot times PER VARIANT -- restrict to a small candidate set.")

  res <- t(sapply(seq_len(ncol(geno.mat)), function(j) {
    g_col <- geno.mat[, j]
    beta_full <- beta_g_curve(seq_len(n), g_col)
    if (any(is.na(beta_full))) return(c(Sn = NA_real_, b = b, pval = NA_real_))
    r_hat <- mean(beta_full)
    v_n <- beta_full - r_hat

    diff_mat <- t(sapply(seq_len(n_boot), function(i) {
      idx <- sample(n, size = b, replace = FALSE)
      beta_sub <- beta_g_curve(idx, g_col)
      beta_sub - beta_full
    }))
    var_tau <- apply(diff_mat, 2, stats::var, na.rm = TRUE)
    var_tau[!is.finite(var_tau) | var_tau <= 0] <- .Machine$double.eps

    if (statistic == "cms") {
      Sn_obs <- (n / b) * mean(v_n^2 / var_tau)
      Sn_boot <- rowMeans(diff_mat^2 / matrix(var_tau, nrow = n_boot, ncol = length(taus), byrow = TRUE))
    } else {
      Sn_obs <- sqrt(n / b) * max(abs(v_n) / sqrt(var_tau))
      Sn_boot <- apply(abs(diff_mat) / matrix(sqrt(var_tau), nrow = n_boot, ncol = length(taus), byrow = TRUE), 1, max)
    }
    pval <- (1 + sum(Sn_boot >= Sn_obs, na.rm = TRUE)) / (1 + sum(!is.na(Sn_boot)))
    c(Sn = Sn_obs, b = b, pval = pval)
  }))

  data.frame(variant_id = variant_ids,
             subsamp_Sn = res[, "Sn"], subsamp_b = res[, "b"], subsamp_pval = res[, "pval"],
             stringsAsFactors = FALSE)
}

#' Quantile TWAS Weight Pipeline
#'
#' @param X Matrix of genotypes
#' @param Y Matrix or vector of phenotypes
#' @param Z Matrix of covariates (optional)
#' @param maf Vector of minor allele frequencies (optional)
#' @param region_id Name of the region being analyzed
#' @param quantile_qtl_tau_list Vector of quantiles for QTL analysis
#' @param quantile_twas_tau_list Vector of quantiles for TWAS analysis
#'
#'
#' @return A list containing various results from the TWAS weight pipeline:
#' \itemize{
#'   \item qr_screen_pvalue_df: Data frame with QR screening results: pavlue, qvalue and zscore.
#'   \item message: Any informational or warning messages.
#'   \item twas_variant_names: Names of variants used in TWAS weight calculation.
#'   \item rq_coef_df: Data frame with quantile regression coefficients.
#'   \item twas_weight: Matrix of TWAS weights.
#'   \item pseudo_R2: Vector of pseudo R-squared values.
#'   \item quantile_twas_prediction: Matrix of TWAS predictions.
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. QR screening to identify significant SNPs.
#' 2. Filtering of highly correlated SNPs.
#' 3. LD clumping and pruning(use filtered SNPs from step 1).
#' 4. Calculation of QR coefficients for selected SNPs(use filtered SNPs from step 3).
#' 5. Calculation of TWAS weights and pseudo R-squared values(use filtered SNPs from step 2).
#'
#' @examples
#' # Example usage:
#' # X <- matrix of genotypes
#' # Y <- vector of phenotypes
#' # Z <- matrix of covariates
#' # results <- quantile_twas_weight_pipeline(X, Y, Z, region_id = "GeneA")
#'
#' @export
quantile_twas_weight_pipeline <- function(X, Y, Z = NULL, maf = NULL, region_id = "",
                                          ld_reference_meta_file = NULL, twas_maf_cutoff = 0.01,
                                          ld_clumping = FALSE, ld_pruning = FALSE,
                                          screen_significant = TRUE,
                                          quantile_qtl_tau_list = seq(0.05, 0.95, by = 0.05),
                                          quantile_twas_tau_list = seq(0.01, 0.99, by = 0.01),
                                          screen_method = "qvalue",
                                          screen_threshold = 0.05,
                                          xi_tau_range = seq(0.1, 0.9, by = 0.05),
                                          heterogeneity_methods = c("wald"),
                                          heterogeneity_n_boot = 199,
                                          keep_variants = NULL,
                                          marginal_beta_calculate = TRUE,
                                          twas_weight_calculate = TRUE,
                                          qrank_screen_calculate = TRUE,
                                          vqtl_calculate = TRUE,
                                          pre_filter_by_pqr = FALSE,
                                          initial_corr_filter_cutoff = 0.8,
                                          full_rank_corr_filter_cutoff = seq(0.75, 0.5, by = -0.05)) {
  # Initialize results list
  results <- list()

  # Step 1: vQTL (optional)
  if (vqtl_calculate) {
    # Step 1-1: Calculate vQTL rank scores
    message("Step 0: Calculating vQTL rank scores for region ", region_id)
    num_tau_levels <- length(quantile_qtl_tau_list) # Convert tau.list to numeric count
    rank_score <- QUAIL_rank_score_pipeline(
      phenotype = Y,
      covariates = Z,
      num_tau_levels = num_tau_levels,
      method = "equal",
      num_cores = 1
    )
    message("vQTL rank scores calculated.")

    # Step 1-2: Run vQTL pipeline
    message("Step 0.5: Running vQTL analysis for rank scores in region ", region_id)
    vqtl_results <- QUAIL_pipeline(
      genotype = X,
      rank_score = rank_score,
      covariates = Z,
      phenotype_id = colnames(Y)[1]
    )
    message("vQTL analysis completed.")
    results$vqtl_results <- vqtl_results
  } else {
    message("Skipping vQTL calculation.")
  }

  if (qrank_screen_calculate) {
    # Step 2: QR screen
    message("Starting QR screen for region ", region_id)
    p.screen <- qr_screen(X = X, Y = Y, Z = Z, tau.list = quantile_qtl_tau_list, screen_threshold = screen_threshold, screen_method = screen_method, top_count = 10, top_percent = 15)
    message(paste0("Number of SNPs after QR screening: ", length(p.screen$sig_SNP_threshold)))
    message("QR screen completed. Screening significant SNPs")
    results$qr_screen_pvalue_df <- p.screen$df_result

    if (screen_significant && length(p.screen$sig_SNP_threshold) == 0) {
      results$message <- paste0("No significant SNPs detected in region ", region_id)
      return(results)
    }

    if (screen_significant) {
      X_filtered <- X[, p.screen$sig_SNP_threshold, drop = FALSE]
    } else {
      X_filtered <- X
    }

    # # Step 3: Optional LD clumping and pruning from results of QR_screen (using original QR screen results)
    if (ld_clumping) {
      message("Performing LD clumping and pruning from QR screen results...")
      LD_SNPs <- multicontext_ld_clumping(X = X[, p.screen$sig_SNP_threshold, drop = FALSE], qr_results = p.screen, maf_list = NULL)
      selected_snps <- if (ld_pruning) LD_SNPs$final_SNPs else LD_SNPs$clumped_SNPs
      x_clumped <- X[, p.screen$sig_SNP_threshold, drop = FALSE][, selected_snps, drop = FALSE]
    } else {
      message("Skipping LD clumping.")
    }

  } else {
    message("Skipping QR screen.")
  }

  # Determine whether to skip marginal beta calculation:
  # - skip if marginal_beta_calculate = FALSE
  # - skip if keep_variants is provided but empty (length 0)
  # - skip if qrank_screen_calculate = FALSE and keep_variants is NULL (no variants to select)
  skip_marginal_beta <- !marginal_beta_calculate ||
    (!is.null(keep_variants) && length(keep_variants) == 0) ||
    (!qrank_screen_calculate && is.null(keep_variants))

  if (!skip_marginal_beta) {
    # Step 4: Fit marginal QR to get beta with SNPs for quantile_qtl_tau_list values
    message("Fitting marginal QR for selected SNPs...")
    if (qrank_screen_calculate) {
      X_for_qr <- if (ld_clumping) x_clumped else X_filtered
      if (!is.null(keep_variants)) {
        variants_to_keep <- intersect(keep_variants, colnames(X_for_qr))
        if (length(variants_to_keep) > 0) {
          X_for_qr <- X_for_qr[, variants_to_keep, drop = FALSE]
          message("Filtered to ", ncol(X_for_qr), " variants from keep_variants list for QR analysis")
        } else {
          message("Warning: No variants from keep_variants found in selected SNPs, using all selected SNPs")
        }
      }
    } else {
      # qrank_screen_calculate = FALSE but keep_variants provided
      variants_to_keep <- intersect(keep_variants, colnames(X))
      if (length(variants_to_keep) > 0) {
        X_for_qr <- X[, variants_to_keep, drop = FALSE]
        message("Using ", ncol(X_for_qr), " variants from keep_variants list for QR analysis (QR screen skipped)")
      } else {
        message("Warning: No variants from keep_variants found in X, skipping marginal beta calculation")
        skip_marginal_beta <- TRUE
      }
    }
  }

  if (!skip_marginal_beta) {
    rq_coef_result <- perform_qr_analysis(X = X_for_qr, Y = Y, Z = Z, tau_values = quantile_qtl_tau_list)

    # Step 5: Heterogeneity calculation
    # Step 5-1: beta_heterogeneity index in marginal model
    message("Marginal QR for selected SNPs completed. Calculating beta heterogeneity...")
    beta_heterogeneity <- calculate_coef_heterogeneity(rq_coef_result)
    message("Beta heterogeneity calculation completed.")

    # Step 5-2: Calculate xi correlation (Chatterjee correlation test)
    message("Calculating xi correlation for QR coefficients...")
    xi_correlation <- calculate_xi_correlation(rq_coef_result, tau_range = xi_tau_range, min_valid = 10)
    message("Xi correlation calculation completed.")

    # Merge xi and xi_pval into rq_coef_result (using left_join to preserve row order)
    rq_coef_result <- rq_coef_result %>%
      dplyr::left_join(xi_correlation, by = "variant_id")

    results$rq_coef_df <- rq_coef_result
    results$beta_heterogeneity <- beta_heterogeneity
    results$xi_correlation <- xi_correlation

    # Step 5-3 (optional): bootstrap-calibrated xi and/or KhmaladzeTest heterogeneity tests.
    # Both need individual-level X_for_qr/Y/Z (not just the saved coef_qr_* columns) and are
    # much more expensive than Step 5-2 -- only requested when heterogeneity_methods asks for
    # them, which should be paired with keep_variants restricted to a small candidate set.
    if ("xi_calibrated" %in% heterogeneity_methods) {
      message("Calculating bootstrap-calibrated xi correlation for QR coefficients...")
      xi_calibrated <- calculate_xi_correlation_calibrated(
        X = X_for_qr, Y = Y, Z = Z, taus = quantile_qtl_tau_list, n_boot = heterogeneity_n_boot
      )
      rq_coef_result <- rq_coef_result %>%
        dplyr::left_join(xi_calibrated, by = "variant_id")
      results$rq_coef_df <- rq_coef_result
      results$xi_calibrated <- xi_calibrated
      message("Bootstrap-calibrated xi correlation completed.")
    }

    if ("khmaladze" %in% heterogeneity_methods) {
      message("Calculating KhmaladzeTest heterogeneity statistics...")
      khmaladze_result <- calculate_khmaladze_heterogeneity(
        X = X_for_qr, Y = Y, Z = Z, taus = quantile_qtl_tau_list, n_boot = heterogeneity_n_boot
      )
      rq_coef_result <- rq_coef_result %>%
        dplyr::left_join(khmaladze_result, by = "variant_id")
      results$rq_coef_df <- rq_coef_result
      results$khmaladze_result <- khmaladze_result
      message("KhmaladzeTest heterogeneity calculation completed.")
    }

    if ("wald" %in% heterogeneity_methods) {
      message("Calculating Wald heterogeneity statistics...")
      wald_result <- calculate_wald_heterogeneity(
        X = X_for_qr, Y = Y, Z = Z, taus = quantile_qtl_tau_list
      )
      rq_coef_result <- rq_coef_result %>%
        dplyr::left_join(wald_result, by = "variant_id")
      results$rq_coef_df <- rq_coef_result
      results$wald_result <- wald_result
      message("Wald heterogeneity calculation completed.")
    }

    if ("subsampling" %in% heterogeneity_methods) {
      message("Calculating subsampling heterogeneity statistics...")
      subsampling_result <- calculate_subsampling_heterogeneity(
        X = X_for_qr, Y = Y, Z = Z, taus = quantile_qtl_tau_list, n_boot = heterogeneity_n_boot
      )
      rq_coef_result <- rq_coef_result %>%
        dplyr::left_join(subsampling_result, by = "variant_id")
      results$rq_coef_df <- rq_coef_result
      results$subsampling_result <- subsampling_result
      message("Subsampling heterogeneity calculation completed.")
    }
  } else {
    message("Skipping marginal beta calculation and heterogeneity analysis.")
  }

  if (twas_weight_calculate && qrank_screen_calculate) {
  # Step 6: Optional LD panel filtering and MAF filtering from results of QR_screen
  if (!is.null(ld_reference_meta_file)) {
    message("Starting LD panel filtering...")
    ld_result <- tryCatch(
      {
        variants_kept <- filter_variants_by_ld_reference(colnames(X_filtered), ld_reference_meta_file)
        if (length(variants_kept$data) == 0) NULL else variants_kept
      },
      error = function(e) {
        message("Error in LD filtering for region ", region_id, ": ", e$message)
        NULL
      }
    )

    if (is.null(ld_result)) {
      results$message <- paste0("No SNPs left or error in LD filtering in region ", region_id)
      return(results)
    }

    X_filtered <- X_filtered[, ld_result$data, drop = FALSE]
    message(paste0("Number of SNPs after LD filtering: ", ncol(X_filtered)))

    # MAF filtering
    if (!is.null(maf)) {
      maf_filtered <- maf[colnames(X_filtered)] > twas_maf_cutoff
      X_filtered <- X_filtered[, maf_filtered, drop = FALSE]

      # Check if any SNPs are left after MAF filtering
      if (ncol(X_filtered) == 0) {
        results$message <- paste0("No SNPs left after MAF filtering in region ", region_id)
        return(results)
      }

      message(paste0("Number of SNPs after MAF filtering: ", ncol(X_filtered)))
    }
  }

    # Step 7: Optionally pre-filter variants by raw p_qr
    if (pre_filter_by_pqr) {
      pqr_attempts <- list(
        list(pval = 0.05, full_rank_corr = full_rank_corr_filter_cutoff),
        list(pval = 0.01, full_rank_corr = full_rank_corr_filter_cutoff)
      )
    } else {
      message("Skipping p_qr pre-filtering, using X_filtered directly...")
      pqr_attempts <- list(
        list(pval = NULL, full_rank_corr = full_rank_corr_filter_cutoff)
      )
    }

    qr_beta_R2_results <- NULL
    for (attempt in seq_along(pqr_attempts)) {
      pqr_cutoff <- pqr_attempts[[attempt]]$pval
      full_rank_corr <- pqr_attempts[[attempt]]$full_rank_corr

      # Pre-filter by p_qr if enabled
      if (!is.null(pqr_cutoff)) {
        raw_pvals <- p.screen$pvec[colnames(X_filtered)]
        sig_raw <- which(raw_pvals < pqr_cutoff)
        if (length(sig_raw) == 0) {
          message("No variants with raw p_qr < ", pqr_cutoff, " in region ", region_id)
          next
        }
        if (length(sig_raw) < ncol(X_filtered)) {
          message("Pre-filtering variants by raw p_qr < ", pqr_cutoff, ": keeping ", length(sig_raw), " out of ", ncol(X_filtered), " variants")
          X_for_corr <- X_filtered[, sig_raw, drop = FALSE]
        } else {
          X_for_corr <- X_filtered
        }
      } else {
        X_for_corr <- X_filtered
      }

      # Step 8: Initial filter of highly correlated SNPs
      message("Filtering highly correlated SNPs...")
      if (ncol(X_for_corr) > 1) {
        filtered <- ld_prune_by_correlation(X_for_corr, initial_corr_filter_cutoff)
        X.filter <- filtered$X.new
      } else {
        X.filter <- X_for_corr
      }

      # Step 9: Fit QR and get twas weight and R2 for all taus
      message("Fitting full QR to calculate TWAS weights and pseudo R-squared values...")
      AssocData <- list(X = X, Y = Y, C = Z, X.filter = X.filter)
      qr_beta_R2_results <- tryCatch(
        {
          calculate_qr_and_pseudo_R2(AssocData, quantile_twas_tau_list, corr_thresholds = full_rank_corr)
        },
        error = function(e) {
          message("Attempt ", attempt, " failed: ", e$message)
          NULL
        }
      )

      if (!is.null(qr_beta_R2_results)) break
    }

    if (is.null(qr_beta_R2_results)) {
      results$message <- paste0("Failed to fit QR model after all attempts in region ", region_id)
      return(results)
    }

  X.filter <- qr_beta_R2_results$X.filter
  message("TWAS weights and pseudo R-squared calculations completed.")

  # Add additional results
  results$twas_variant_names <- colnames(X.filter)
  results$twas_weight <- qr_beta_R2_results$beta_mat
  results$pseudo_R2 <- qr_beta_R2_results$pseudo_R2
  results$quantile_twas_prediction <- X.filter %*% results$twas_weight
  } else {
    message("Skipping TWAS weight calculation.")
  }

  return(results)
}
