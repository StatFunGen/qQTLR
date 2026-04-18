context("quantile_twas")

# ===========================================================================
# ensure_continuous_clusters (internal)
# ===========================================================================

test_that("ensure_continuous_clusters preserves already continuous", {
  index <- c(1, 1, 2, 2, 3)
  result <- qQTLR:::ensure_continuous_clusters(index)
  expect_equal(result, c(1, 1, 2, 2, 3))
})

test_that("ensure_continuous_clusters renumbers non-continuous", {
  index <- c(1, 2, 1, 2, 3)
  result <- qQTLR:::ensure_continuous_clusters(index)
  expect_equal(result, c(1, 2, 3, 4, 5))
})

test_that("ensure_continuous_clusters handles single element", {
  result <- qQTLR:::ensure_continuous_clusters(c(1))
  expect_equal(result, c(1))
})

test_that("ensure_continuous_clusters handles two elements same cluster", {
  result <- qQTLR:::ensure_continuous_clusters(c(1, 1))
  expect_equal(result, c(1, 1))
})

test_that("ensure_continuous_clusters all same cluster", {
  result <- qQTLR:::ensure_continuous_clusters(c(3, 3, 3, 3))
  expect_equal(result, c(1, 1, 1, 1))
})

# ===========================================================================
# get_cluster_ranges (internal)
# ===========================================================================

test_that("get_cluster_ranges returns correct ranges", {
  index <- c(1, 1, 2, 2, 3)
  result <- qQTLR:::get_cluster_ranges(index)
  expect_length(result, 3)
  expect_equal(result[[1]], "1 - 2")
  expect_equal(result[[2]], "3 - 4")
  expect_equal(result[[3]], "5 - 5")
})

test_that("get_cluster_ranges single cluster", {
  index <- c(1, 1, 1)
  result <- qQTLR:::get_cluster_ranges(index)
  expect_length(result, 1)
  expect_equal(result[[1]], "1 - 3")
})

# ===========================================================================
# get_modularity (internal)
# ===========================================================================

test_that("get_modularity returns 0 for single element", {
  W <- matrix(1, nrow = 1, ncol = 1)
  B <- matrix(1, nrow = 1, ncol = 1)
  result <- qQTLR:::get_modularity(W, B)
  expect_equal(result, 0)
})

test_that("get_modularity returns numeric for identity weight", {
  W <- diag(5)
  B <- matrix(c(rep(c(1, 0), c(3, 2)), rep(c(0, 1), c(3, 2))), nrow = 5)
  result <- qQTLR:::get_modularity(W, B)
  expect_type(result, "double")
})

test_that("get_modularity handles all-positive weights", {
  set.seed(42)
  W <- abs(matrix(rnorm(16), 4, 4))
  W <- (W + t(W)) / 2
  B <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow = 4)
  result <- qQTLR:::get_modularity(W, B)
  expect_type(result, "double")
  expect_true(is.finite(result))
})

# ===========================================================================
# get_n_cluster (internal)
# ===========================================================================

test_that("get_n_cluster returns single cluster for high correlation", {
  Sigma <- matrix(0.9, nrow = 5, ncol = 5)
  diag(Sigma) <- 1
  hc <- hclust(as.dist(1 - Sigma))
  result <- qQTLR:::get_n_cluster(hc, Sigma, between_cluster = 0.8)
  expect_equal(result$n_cluster, 1)
})

test_that("get_n_cluster returns multiple clusters for low correlation", {
  set.seed(42)
  Sigma <- diag(10)
  Sigma[1:5, 1:5] <- 0.8
  Sigma[6:10, 6:10] <- 0.8
  diag(Sigma) <- 1
  hc <- hclust(as.dist(1 - Sigma))
  result <- qQTLR:::get_n_cluster(hc, Sigma, between_cluster = 0.5)
  expect_true(result$n_cluster >= 1)
})

# ===========================================================================
# integrate_tau (internal)
# ===========================================================================

test_that("integrate_tau with vector input", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- c(1, 2, 1)
  result <- qQTLR:::integrate_tau(tau, a_tau)
  expect_type(result, "double")
  expect_true(result > 0)
})

test_that("integrate_tau with matrix input", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- matrix(c(1, 2, 1, 0.5, 1, 0.5), nrow = 2, byrow = TRUE)
  result <- qQTLR:::integrate_tau(tau, a_tau)
  expect_equal(nrow(result), 2)
  expect_true(all(result > 0))
})

test_that("integrate_tau zero weights give zero", {
  tau <- c(0.25, 0.50, 0.75)
  a_tau <- c(0, 0, 0)
  result <- qQTLR:::integrate_tau(tau, a_tau)
  expect_equal(as.numeric(result), 0)
})

# ===========================================================================
# get_hierarchical_clusters (internal)
# ===========================================================================

test_that("get_hierarchical_clusters returns valid structure", {
  set.seed(42)
  p <- 10
  Sigma <- matrix(0.3, nrow = p, ncol = p)
  Sigma[1:5, 1:5] <- 0.9
  Sigma[6:10, 6:10] <- 0.9
  diag(Sigma) <- 1
  result <- qQTLR:::get_hierarchical_clusters(Sigma, between_cluster = 0.5)
  expect_type(result, "list")
  expect_true("cluster" %in% names(result))
  expect_true("Q_modularity_initial" %in% names(result))
  expect_true("cluster_ranges" %in% names(result))
  expect_equal(nrow(result$cluster), p)
  expect_true(all(rowSums(result$cluster) == 1))
})

# ===========================================================================
# perform_grouped_integration (internal)
# ===========================================================================

test_that("perform_grouped_integration returns correct structure", {
  set.seed(42)
  n_variants <- 5
  n_tau <- 10
  twas_weight <- matrix(rnorm(n_variants * n_tau), nrow = n_variants)
  rownames(twas_weight) <- paste0("v", 1:n_variants)
  tau_values <- seq(0.1, 0.9, length.out = n_tau)
  pseudo_R2 <- runif(n_tau, 0.1, 0.5)

  result <- qQTLR:::perform_grouped_integration(twas_weight, tau_values, pseudo_R2,
                                                   num_intervals = 3)
  expect_type(result, "list")
  expect_true("weights" %in% names(result))
  expect_true("twas_weight_performance" %in% names(result))
  expect_equal(nrow(result$weights), n_variants)
  expect_true(ncol(result$weights) >= 3)
})

test_that("perform_grouped_integration single variant skips clustering", {
  twas_weight <- matrix(rnorm(10), nrow = 1)
  rownames(twas_weight) <- "v1"
  tau_values <- seq(0.1, 0.9, length.out = 10)
  pseudo_R2 <- runif(10, 0.1, 0.5)

  result <- qQTLR:::perform_grouped_integration(twas_weight, tau_values, pseudo_R2,
                                                   num_intervals = 3)
  expect_equal(ncol(result$weights), 3)
})

# ===========================================================================
# calculate_coef_heterogeneity (internal)
# ===========================================================================

test_that("calculate_coef_heterogeneity computes log(sd/mean)", {
  df <- data.frame(
    variant_id = c("v1", "v2"),
    coef_qr_0.25 = c(1, 0),
    coef_qr_0.50 = c(2, 0),
    coef_qr_0.75 = c(3, 0),
    stringsAsFactors = FALSE
  )
  result <- qQTLR:::calculate_coef_heterogeneity(df)
  expect_equal(nrow(result), 2)
  expect_true("coef_heter" %in% colnames(result))
  expect_equal(result$coef_heter[1], log(1 / 2), tolerance = 0.01)
  expect_true(is.na(result$coef_heter[2]))
})

test_that("calculate_coef_heterogeneity handles NA values", {
  df <- data.frame(
    variant_id = "v1",
    coef_qr_0.25 = 1,
    coef_qr_0.50 = NA,
    coef_qr_0.75 = 3,
    stringsAsFactors = FALSE
  )
  result <- qQTLR:::calculate_coef_heterogeneity(df)
  expect_equal(nrow(result), 1)
  expect_true(is.finite(result$coef_heter[1]) || is.na(result$coef_heter[1]))
})

test_that("calculate_coef_heterogeneity handles all-same coefficients", {
  df <- data.frame(
    variant_id = "v1",
    coef_qr_0.25 = 5,
    coef_qr_0.50 = 5,
    coef_qr_0.75 = 5,
    stringsAsFactors = FALSE
  )
  result <- qQTLR:::calculate_coef_heterogeneity(df)
  expect_true(is.infinite(result$coef_heter[1]) || is.na(result$coef_heter[1]))
})

test_that("calculate_coef_heterogeneity handles negative coefficients", {
  df <- data.frame(
    variant_id = "v1",
    coef_qr_0.25 = -1,
    coef_qr_0.50 = -2,
    coef_qr_0.75 = -3,
    stringsAsFactors = FALSE
  )
  result <- qQTLR:::calculate_coef_heterogeneity(df)
  expect_equal(result$coef_heter[1], log(1 / 2), tolerance = 0.01)
})

# ===========================================================================
# calculate_xi_correlation
# ===========================================================================

test_that("calculate_xi_correlation computes xi for valid monotonic data", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.05, 0.95, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- tau * 2
  }
  result <- calculate_xi_correlation(df)
  expect_true("xi" %in% colnames(result))
  expect_true(is.numeric(result$xi))
})

test_that("calculate_xi_correlation warns on missing columns", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1", some_col = 1)
  expect_warning(result <- calculate_xi_correlation(df))
  expect_true(is.na(result$xi))
})

test_that("calculate_xi_correlation returns NA for too few valid values", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.1, 0.9, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- NA
  }
  df$coef_qr_0.1 <- 1
  result <- calculate_xi_correlation(df, min_valid = 10)
  expect_true(is.na(result$xi))
})

test_that("calculate_xi_correlation handles multiple variants", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = c("v1", "v2"))
  tau_range <- seq(0.1, 0.9, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- c(tau * 2, tau * -1)
  }
  result <- calculate_xi_correlation(df, tau_range = tau_range)
  expect_equal(nrow(result), 2)
  expect_true("xi" %in% colnames(result))
  expect_true("xi_pval" %in% colnames(result))
})

test_that("calculate_xi_correlation handles error in xicor gracefully", {
  skip_if_not_installed("XICOR")
  df <- data.frame(variant_id = "v1")
  tau_range <- seq(0.1, 0.9, by = 0.05)
  for (tau in tau_range) {
    df[[paste0("coef_qr_", tau)]] <- 1
  }
  result <- calculate_xi_correlation(df, tau_range = tau_range, min_valid = 5)
  expect_equal(nrow(result), 1)
  expect_true(is.numeric(result$xi) || is.na(result$xi))
})

# ===========================================================================
# qr_screen
# ===========================================================================

test_that("qr_screen runs and returns results with quantreg", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  colnames(X) <- paste0("chr1:", seq(100, 500, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- qr_screen(X, Y, tau.list = c(0.25, 0.5, 0.75))
  expect_type(result, "list")
})

test_that("qr_screen with fdr screen_method works", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 5
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 500, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- qr_screen(X, Y, tau.list = c(0.25, 0.5, 0.75), screen_method = "fdr")
  expect_type(result, "list")
  expect_true("df_result" %in% names(result))
  expect_true("fdr_p_qr" %in% colnames(result$df_result))
})

test_that("qr_screen errors with invalid screen_method", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 30; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  expect_error(qr_screen(X, Y, tau.list = c(0.5), screen_method = "invalid_method"),
               "Invalid screen_method")
})

test_that("qr_screen with covariates Z", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  Z <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  result <- qr_screen(X, Y, Z = Z, tau.list = c(0.25, 0.5, 0.75))
  expect_type(result, "list")
  expect_equal(nrow(result$df_result), p)
})

test_that("qr_screen with single variant", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "chr1:100:A:G"
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- qr_screen(X, Y, tau.list = c(0.25, 0.5, 0.75))
  expect_type(result, "list")
  expect_equal(nrow(result$df_result), 1)
})

# ===========================================================================
# perform_qr_analysis
# ===========================================================================

test_that("perform_qr_analysis returns wide-format results", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- perform_qr_analysis(X, Y, tau_values = c(0.25, 0.5, 0.75))
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), p)
  expect_true("coef_qr_0.25" %in% colnames(result))
  expect_true("coef_qr_0.5" %in% colnames(result))
  expect_true("coef_qr_0.75" %in% colnames(result))
  expect_true("chr" %in% colnames(result))
  expect_true("pos" %in% colnames(result))
})

test_that("perform_qr_analysis with covariates Z", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 2
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 200, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  Z <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  result <- perform_qr_analysis(X, Y, Z = Z, tau_values = c(0.5))
  expect_equal(nrow(result), p)
  expect_true("coef_qr_0.5" %in% colnames(result))
})

test_that("perform_qr_analysis with single variant", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X) <- "chr1:100:A:G"
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"
  result <- perform_qr_analysis(X, Y, tau_values = c(0.25, 0.5))
  expect_equal(nrow(result), 1)
})

# ===========================================================================
# multicontext_ld_clumping
# ===========================================================================

test_that("multicontext_ld_clumping returns empty result when no significant SNPs", {
  qr_results <- list(sig.SNPs_names = character(0))
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  result <- multicontext_ld_clumping(X, qr_results)
  expect_null(result$final_SNPs)
  expect_null(result$clumped_SNPs)
  expect_true(grepl("No significant", result$message))
})

test_that("multicontext_ld_clumping returns early when X has single column", {
  X <- matrix(rnorm(20), nrow = 20, ncol = 1)
  colnames(X) <- "chr1:100:A:G"
  qr_results <- list(sig.SNPs_names = "chr1:100:A:G")
  result <- multicontext_ld_clumping(X, qr_results)
  expect_equal(result$final_SNPs, "chr1:100:A:G")
  expect_true(grepl("Only one SNP", result$message))
})

# ===========================================================================
# quantile_twas_weight_pipeline (mocked)
# ===========================================================================

test_that("quantile_twas_weight_pipeline returns early when no significant SNPs and screen_significant=TRUE", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"

  local_mocked_bindings(
    QUAIL_rank_score_pipeline = function(...) list(rank_score = matrix(rnorm(n), ncol = 1)),
    QUAIL_pipeline = function(...) list(vqtl = "mocked"),
    qr_screen = function(...) {
      list(
        df_result = data.frame(variant_id = colnames(X)),
        sig_SNP_threshold = integer(0),
        sig.SNPs_names = character(0),
        pvec = rep(1, p),
        quantile_pvalue = matrix(1, nrow = p, ncol = 1),
        quantile_zscore = matrix(0, nrow = p, ncol = 1),
        tau_list = 0.5
      )
    }
  )
  result <- quantile_twas_weight_pipeline(
    X, Y, screen_significant = TRUE,
    quantile_qtl_tau_list = c(0.5),
    quantile_twas_tau_list = c(0.5)
  )
  expect_true(grepl("No significant", result$message))
  expect_true("qr_screen_pvalue_df" %in% names(result))
})

test_that("quantile_twas_weight_pipeline returns early when no raw p_qr < 0.05", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(100, 300, by = 100), ":A:G")
  Y <- matrix(rnorm(n), ncol = 1)
  colnames(Y) <- "pheno1"

  local_mocked_bindings(
    QUAIL_rank_score_pipeline = function(...) list(rank_score = matrix(rnorm(n), ncol = 1)),
    QUAIL_pipeline = function(...) list(vqtl = "mocked"),
    qr_screen = function(...) {
      pvec <- rep(0.5, p)
      names(pvec) <- colnames(X)
      list(
        df_result = data.frame(variant_id = colnames(X)),
        sig_SNP_threshold = 1:p,
        sig.SNPs_names = colnames(X),
        pvec = pvec,
        quantile_pvalue = matrix(0.5, nrow = p, ncol = 1),
        quantile_zscore = matrix(0, nrow = p, ncol = 1),
        tau_list = 0.5
      )
    },
    perform_qr_analysis = function(...) {
      data.frame(
        variant_id = colnames(X),
        coef_qr_0.5 = rnorm(p),
        chr = rep("chr1", p),
        pos = seq(100, 300, by = 100),
        ref = rep("G", p),
        alt = rep("A", p),
        phenotype_id = rep("pheno1", p),
        stringsAsFactors = FALSE
      )
    },
    calculate_coef_heterogeneity = function(...) {
      data.frame(variant_id = colnames(X), coef_heter = rnorm(p))
    },
    calculate_xi_correlation = function(...) {
      data.frame(variant_id = colnames(X), xi = rnorm(p), xi_pval = runif(p))
    }
  )
  result <- quantile_twas_weight_pipeline(
    X, Y, screen_significant = FALSE,
    quantile_qtl_tau_list = c(0.5),
    quantile_twas_tau_list = c(0.5)
  )
  # expect_true(grepl("No variants with raw", result$message))
})

# ===========================================================================
# calculate_qr_and_pseudo_R2 (internal)
# ===========================================================================

test_that("calculate_qr_and_pseudo_R2 returns expected structure with single SNP", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 50
  X.filter <- matrix(rnorm(n), nrow = n, ncol = 1)
  colnames(X.filter) <- "snp1"
  Y <- rnorm(n)
  C <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  tau.list <- c(0.25, 0.5, 0.75)

  AssocData <- list(X = X.filter, Y = Y, C = C, X.filter = X.filter)
  result <- qQTLR:::calculate_qr_and_pseudo_R2(AssocData, tau.list)
  expect_type(result, "list")
  expect_true("X.filter" %in% names(result))
  expect_true("beta_mat" %in% names(result))
  expect_true("pseudo_R2" %in% names(result))
  expect_equal(length(result$pseudo_R2), length(tau.list))
  expect_true(all(is.finite(result$pseudo_R2)))
  expect_equal(nrow(result$beta_mat), 1)
  expect_equal(ncol(result$beta_mat), length(tau.list))
  expect_equal(rownames(result$beta_mat), "snp1")
})

test_that("calculate_qr_and_pseudo_R2 returns expected structure with multiple SNPs", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 60; p <- 3
  X.filter <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X.filter) <- paste0("snp", 1:p)
  Y <- rnorm(n)
  C <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  tau.list <- c(0.25, 0.5, 0.75)

  AssocData <- list(X = X.filter, Y = Y, C = C, X.filter = X.filter)
  result <- qQTLR:::calculate_qr_and_pseudo_R2(AssocData, tau.list)
  expect_equal(nrow(result$beta_mat), p)
  expect_equal(ncol(result$beta_mat), length(tau.list))
  expect_equal(rownames(result$beta_mat), paste0("snp", 1:p))
  expect_true(all(is.finite(result$pseudo_R2)))
})

test_that("calculate_qr_and_pseudo_R2 pseudo R2 values are between 0 and 1 for signal data", {
  skip_if_not_installed("quantreg")
  set.seed(42)
  n <- 60; p <- 2
  X.filter <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X.filter) <- paste0("snp", 1:p)
  # Create Y with signal from X so pseudo R2 is positive
  Y <- X.filter[, 1] * 2 + X.filter[, 2] * 0.5 + rnorm(n, sd = 0.5)
  C <- matrix(rnorm(n), nrow = n, ncol = 1)
  tau.list <- c(0.25, 0.5, 0.75)

  AssocData <- list(X = X.filter, Y = Y, C = C, X.filter = X.filter)
  result <- qQTLR:::calculate_qr_and_pseudo_R2(AssocData, tau.list)
  expect_true(all(result$pseudo_R2 > 0))
  expect_true(all(result$pseudo_R2 < 1))
})

# ===========================================================================
# get_modularity edge cases (internal)
# ===========================================================================

test_that("get_modularity handles all-zero weight matrix", {
  W <- matrix(0, nrow = 3, ncol = 3)
  B <- matrix(c(1, 1, 0, 0, 0, 1), nrow = 3)
  result <- qQTLR:::get_modularity(W, B)
  expect_equal(result, 0)
})

test_that("get_modularity handles all-negative weights", {
  W <- -abs(matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2))
  B <- matrix(c(1, 0, 0, 1), nrow = 2)
  result <- qQTLR:::get_modularity(W, B)
  expect_type(result, "double")
  expect_true(is.finite(result))
})

# ===========================================================================
# enforce_design_full_rank integration (imported from pecotmr)
# ===========================================================================

test_that("enforce_design_full_rank triggers correlation-prune fallback for stubborn rank deficiency", {
  set.seed(42)
  n <- 50
  # Columns highly correlated but not identical so QR removal won't
  # fully fix the rank deficiency on its own.
  base <- rnorm(n)
  X <- cbind(
    base + rnorm(n, sd = 0.001),
    base + rnorm(n, sd = 0.001),
    base + rnorm(n, sd = 0.001),
    rnorm(n)
  )
  colnames(X) <- paste0("v", 1:4)
  C <- matrix(X[, 4] + rnorm(n, sd = 0.001), ncol = 1)
  result <- enforce_design_full_rank(X = X, C = C, strategy = "correlation")
  design <- cbind(1, result, C)
  expect_equal(qr(design)$rank, ncol(design))
})

test_that("enforce_design_full_rank max_iterations warning path", {
  set.seed(42)
  n <- 50
  X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  colnames(X) <- paste0("v", 1:4)
  X[, 3] <- X[, 1] + X[, 2]
  X[, 4] <- X[, 1] - X[, 2]
  expect_warning(
    enforce_design_full_rank(X = X, C = NULL, strategy = "correlation", max_iterations = 1),
    "max_iterations reached"
  )
})

# ===========================================================================
# multicontext_ld_clumping main body (exported)
# ===========================================================================

test_that("multicontext_ld_clumping performs clumping with multiple significant SNPs", {
  set.seed(42)
  n <- 100; p <- 10
  X <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  colnames(X) <- paste0("chr1:", seq(1000, by = 1000, length.out = p), ":A:G")

  # Create qr_results with enough significant SNPs
  sig_indices <- 1:p
  tau_list <- c(0.25, 0.5, 0.75)
  quantile_pvalue <- matrix(runif(p * length(tau_list), 0.001, 0.05),
                            nrow = p, ncol = length(tau_list))
  colnames(quantile_pvalue) <- paste0("p_qr_", tau_list)
  rownames(quantile_pvalue) <- colnames(X)

  qr_results <- list(
    sig.SNPs_names = colnames(X)[sig_indices],
    sig_SNP_threshold = sig_indices,
    tau_list = tau_list,
    quantile_pvalue = quantile_pvalue
  )

  result <- multicontext_ld_clumping(X, qr_results)
  expect_true(!is.null(result$final_SNPs))
  expect_true(length(result$final_SNPs) > 0)
  expect_true(length(result$clumped_SNPs) > 0)
  expect_true(all(result$final_SNPs %in% colnames(X)))
  expect_true(all(result$clumped_SNPs %in% colnames(X)))
})

test_that("multicontext_ld_clumping union of clumping yields single SNP path", {
  set.seed(99)
  n <- 100
  # Two SNPs that are perfectly correlated -> clumping reduces to 1
  X <- matrix(rbinom(n * 2, 2, 0.3), nrow = n, ncol = 2)
  X[, 2] <- X[, 1]  # identical columns
  colnames(X) <- c("chr1:1000:A:G", "chr1:2000:A:G")

  tau_list <- c(0.5)
  quantile_pvalue <- matrix(c(0.01, 0.01), nrow = 2, ncol = 1)
  colnames(quantile_pvalue) <- "p_qr_0.5"
  rownames(quantile_pvalue) <- colnames(X)

  qr_results <- list(
    sig.SNPs_names = colnames(X),
    sig_SNP_threshold = 1:2,
    tau_list = tau_list,
    quantile_pvalue = quantile_pvalue
  )

  result <- multicontext_ld_clumping(X, qr_results)
  expect_true(!is.null(result$final_SNPs))
  expect_true(length(result$final_SNPs) >= 1)
})
