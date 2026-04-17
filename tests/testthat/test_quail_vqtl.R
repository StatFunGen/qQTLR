context("quail_vqtl")

# ===========================================================================
# Helpers: small synthetic genotype / phenotype data
# ===========================================================================

make_genotype_matrix <- function(n = 50, p = 5, seed = 123) {
  set.seed(seed)
  geno <- matrix(rbinom(n * p, 2, 0.3), nrow = n, ncol = p)
  colnames(geno) <- paste0("chr1:", seq(100, by = 100, length.out = p), ":A:G")
  geno
}

make_rank_score <- function(n = 50, seed = 123) {
  set.seed(seed)
  rnorm(n)
}

make_covariates <- function(n = 50, seed = 123) {
  set.seed(seed + 1)
  data.frame(
    age = rnorm(n, 50, 10),
    sex = rbinom(n, 1, 0.5)
  )
}

# ===========================================================================
# QUAIL_pipeline tests
# ===========================================================================

# ---------- Input validation -----------------------------------------------

test_that("QUAIL_pipeline errors when rank_score is not numeric", {
  geno <- make_genotype_matrix()
  expect_error(
    QUAIL_pipeline(geno, rank_score = c("a", "b", "c")),
    "rank_score must be a numeric vector"
  )
})

# ---------- Basic synthetic data run ---------------------------------------

test_that("QUAIL_pipeline runs with basic synthetic data", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  result <- QUAIL_pipeline(geno, rs)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), ncol(geno))
})

# ---------- Output structure -----------------------------------------------

test_that("QUAIL_pipeline output has required columns", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  result <- QUAIL_pipeline(geno, rs)
  expected_cols <- c("phenotype_id", "chr", "pos", "A1", "A2",
                     "variant_id", "beta", "se", "z", "p", "q", "N")
  expect_true(all(expected_cols %in% colnames(result)))
})

test_that("QUAIL_pipeline N column equals number of samples", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix(n = 40)
  rs <- make_rank_score(n = 40)
  result <- QUAIL_pipeline(geno, rs)
  expect_true(all(result$N == 40))
})

test_that("QUAIL_pipeline variant_id matches genotype column names", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  result <- QUAIL_pipeline(geno, rs)
  expect_equal(result$variant_id, colnames(geno))
})

test_that("QUAIL_pipeline parses chr and pos from variant IDs", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  result <- QUAIL_pipeline(geno, rs)
  # All variants are on chr1
  expect_true(all(result$chr == 1))
  # Positions should be 100, 200, ..., 500
  expect_equal(result$pos, seq(100, by = 100, length.out = ncol(geno)))
})

test_that("QUAIL_pipeline phenotype_id defaults to NA when not provided", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  result <- QUAIL_pipeline(geno, rs)
  expect_true(all(is.na(result$phenotype_id)))
})

test_that("QUAIL_pipeline phenotype_id is set when provided", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  result <- QUAIL_pipeline(geno, rs, phenotype_id = "my_trait")
  expect_true(all(result$phenotype_id == "my_trait"))
})

# ---------- With covariates ------------------------------------------------

test_that("QUAIL_pipeline works with covariates", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  cov <- make_covariates()
  result <- QUAIL_pipeline(geno, rs, covariates = cov)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), ncol(geno))
})

test_that("QUAIL_pipeline covariates affect beta estimates", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  cov <- make_covariates()
  result_no_cov <- QUAIL_pipeline(geno, rs)
  result_with_cov <- QUAIL_pipeline(geno, rs, covariates = cov)
  # Betas should differ when adjusting for covariates
  expect_false(identical(result_no_cov$beta, result_with_cov$beta))
})

test_that("QUAIL_pipeline accepts matrix covariates", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  cov_mat <- as.matrix(make_covariates())
  result <- QUAIL_pipeline(geno, rs, covariates = cov_mat)
  expect_true(is.data.frame(result))
})

# ---------- p-values and z-scores ------------------------------------------

test_that("QUAIL_pipeline p-values are between 0 and 1", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix()
  rs <- make_rank_score()
  result <- QUAIL_pipeline(geno, rs)
  expect_true(all(result$p >= 0 & result$p <= 1))
})

test_that("QUAIL_pipeline z-scores are finite", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix(n = 60, p = 4)
  rs <- make_rank_score(n = 60)
  result <- QUAIL_pipeline(geno, rs)
  # z-scores should all be finite for well-conditioned data
  expect_true(all(is.finite(result$z)))
})

# ===========================================================================
# univariate_regression tests (internal function)
# ===========================================================================

test_that("univariate_regression returns expected list structure", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(200), 50, 4)
  colnames(X) <- paste0("snp", 1:4)
  y <- rnorm(50)
  res <- qQTLR:::univariate_regression(X, y)
  expect_true(is.list(res))
  expected_names <- c("betahat", "sebetahat", "z_scores", "p_values", "q_values")
  expect_true(all(expected_names %in% names(res)))
  # Each element should have length equal to ncol(X)
  for (nm in expected_names) {
    expect_equal(length(res[[nm]]), ncol(X))
  }
})

test_that("univariate_regression names outputs by column names of X", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(150), 50, 3)
  colnames(X) <- c("var_A", "var_B", "var_C")
  y <- rnorm(50)
  res <- qQTLR:::univariate_regression(X, y)
  expect_equal(names(res$betahat), colnames(X))
  expect_equal(names(res$sebetahat), colnames(X))
})

test_that("univariate_regression handles centering", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(200, mean = 10), 50, 4)
  colnames(X) <- paste0("s", 1:4)
  y <- rnorm(50, mean = 5)
  # With centering (default)
  res_centered <- qQTLR:::univariate_regression(X, y, center = TRUE)
  # Without centering
  res_no_center <- qQTLR:::univariate_regression(X, y, center = FALSE)
  # Since the regression includes an intercept, centering X and y
  # should yield numerically equivalent slope estimates
  expect_equal(res_centered$betahat, res_no_center$betahat, tolerance = 1e-10)
  expect_equal(res_centered$sebetahat, res_no_center$sebetahat, tolerance = 1e-10)
  # Both results should have valid finite values
  expect_true(all(is.finite(res_centered$betahat)))
  expect_true(all(is.finite(res_no_center$betahat)))
})

test_that("univariate_regression handles scaling", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(200, sd = 5), 50, 4)
  colnames(X) <- paste0("s", 1:4)
  y <- rnorm(50)
  res_scaled <- qQTLR:::univariate_regression(X, y, scale = TRUE)
  res_not_scaled <- qQTLR:::univariate_regression(X, y, scale = FALSE)
  # Scaling divides each column by its SD, so beta_scaled = beta_unscaled * column_sd
  col_sds <- apply(X, 2, sd)
  expect_equal(res_scaled$betahat, res_not_scaled$betahat * col_sds,
               tolerance = 1e-10)
  # Both should be finite
  expect_true(all(is.finite(res_scaled$betahat)))
  expect_true(all(is.finite(res_not_scaled$betahat)))
})

test_that("univariate_regression handles covariates", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(200), 50, 4)
  colnames(X) <- paste0("s", 1:4)
  Z <- matrix(rnorm(100), 50, 2)
  y <- rnorm(50)
  res_cov <- qQTLR:::univariate_regression(X, y, Z = Z)
  res_no_cov <- qQTLR:::univariate_regression(X, y)
  # Adjusting for covariates should generally reduce residual SE
  # (or at minimum keep it comparable)
  expect_true(all(res_cov$sebetahat <= res_no_cov$sebetahat + 1e-10))
  # Both should return valid finite results
  expect_true(all(is.finite(res_cov$betahat)))
  expect_true(all(is.finite(res_cov$sebetahat)))
})

test_that("univariate_regression handles NAs in y", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(200), 50, 4)
  colnames(X) <- paste0("s", 1:4)
  y <- rnorm(50)
  y[c(1, 5, 10)] <- NA
  res <- qQTLR:::univariate_regression(X, y)
  # Should still return valid results
  expect_equal(length(res$betahat), 4)
  expect_true(all(is.finite(res$betahat)))
})

test_that("univariate_regression return_residuals works with covariates", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(200), 50, 4)
  colnames(X) <- paste0("s", 1:4)
  Z <- matrix(rnorm(100), 50, 2)
  y <- rnorm(50)
  res <- qQTLR:::univariate_regression(X, y, Z = Z, return_residuals = TRUE)
  expect_true("residuals" %in% names(res))
  expect_equal(length(res$residuals), length(y))
})

test_that("univariate_regression return_residuals without covariates omits residuals", {
  skip_if_not_installed("qvalue")
  set.seed(42)
  X <- matrix(rnorm(200), 50, 4)
  colnames(X) <- paste0("s", 1:4)
  y <- rnorm(50)
  res <- qQTLR:::univariate_regression(X, y, return_residuals = TRUE)
  # Without Z, residuals should not be in the output
  expect_false("residuals" %in% names(res))
})

# ===========================================================================
# run_linear_regression tests (internal function)
# ===========================================================================

test_that("run_linear_regression returns a data frame with expected columns", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix(n = 40, p = 3)
  pheno <- make_rank_score(n = 40)
  res <- qQTLR:::run_linear_regression(geno, pheno)
  expected_cols <- c("phenotype_id", "chr", "pos", "A1", "A2",
                     "variant_id", "beta", "se", "z", "p", "q", "N")
  expect_true(all(expected_cols %in% colnames(res)))
  expect_equal(nrow(res), ncol(geno))
})

test_that("run_linear_regression passes phenotype_id through", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix(n = 40, p = 3)
  pheno <- make_rank_score(n = 40)
  res <- qQTLR:::run_linear_regression(geno, pheno, phenotype_id = "test_pheno")
  expect_true(all(res$phenotype_id == "test_pheno"))
})

# ===========================================================================
# Edge cases
# ===========================================================================

test_that("QUAIL_pipeline works with a single SNP", {
  skip_if_not_installed("qvalue")
  set.seed(99)
  n <- 30
  geno <- matrix(rbinom(n, 2, 0.3), nrow = n, ncol = 1)
  colnames(geno) <- "chr1:500:C:T"
  rs <- rnorm(n)
  result <- QUAIL_pipeline(geno, rs)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 1)
})

test_that("univariate_regression handles all-zero genotype column gracefully", {
  skip_if_not_installed("qvalue")
  set.seed(99)
  n <- 40
  X <- matrix(rnorm(n * 3), n, 3)
  X[, 2] <- 0 # all-zero column
  colnames(X) <- paste0("s", 1:3)
  y <- rnorm(n)
  # Should not error; the zero-variance column gets NaN which is set to 0
  res <- qQTLR:::univariate_regression(X, y)
  expect_equal(length(res$betahat), 3)
  # The all-zero column (after centering, NaN set to 0) produces beta = 0
  expect_equal(unname(res$betahat[2]), 0)
})

test_that("QUAIL_pipeline se values are non-negative", {
  skip_if_not_installed("qvalue")
  geno <- make_genotype_matrix(n = 60, p = 4)
  rs <- make_rank_score(n = 60)
  result <- QUAIL_pipeline(geno, rs)
  expect_true(all(result$se >= 0))
})
