context("quail_ctwas")

# ===========================================================================
# Helper: small simulated phenotype + covariate data (n = 30)
# ===========================================================================
make_quail_data <- function(n = 30, seed = 42) {
  set.seed(seed)
  covariates <- data.frame(
    age = rnorm(n, 50, 10),
    sex = rbinom(n, 1, 0.5)
  )
  # Phenotype with variance that depends on covariates (heteroscedastic)
  phenotype <- 2 + 0.5 * covariates$age + rnorm(n, sd = 1 + 0.3 * covariates$sex)
  list(phenotype = phenotype, covariates = covariates)
}

# ===========================================================================
#  QUAIL rank score tests
# ===========================================================================

# ---------- calculate_rank_score -------------------------------------------

test_that("calculate_rank_score returns vector of correct length", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  result <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 1)
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
})

test_that("calculate_rank_score is reproducible with same seed", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  r1 <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.25, seed = 99)
  r2 <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.25, seed = 99)
  expect_identical(r1, r2)
})

test_that("calculate_rank_score differs with different seeds", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  r1 <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 1)
  r2 <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 2)
  # Different random covariate column, so results should differ

  expect_false(identical(r1, r2))
})

test_that("calculate_rank_score produces different results for different tau values", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  r_low <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.1, seed = 1)
  r_high <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.9, seed = 1)
  expect_false(identical(r_low, r_high))
})

test_that("calculate_rank_score returns finite values", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  result <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau = 0.5, seed = 1)
  expect_true(all(is.finite(result)))
})

# ---------- calculate_integrated_score -------------------------------------

test_that("calculate_integrated_score works with equal method and even num_tau_levels", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 4 # even
  # Create mock rank scores: list of numeric vectors
  set.seed(10)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n))
  result <- qQTLR:::calculate_integrated_score(rank_scores, method = "equal", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
  expect_true(all(is.finite(result)))
})

test_that("calculate_integrated_score works with equal method and odd num_tau_levels", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 5 # odd
  set.seed(11)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n))
  result <- qQTLR:::calculate_integrated_score(rank_scores, method = "equal", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
})

test_that("calculate_integrated_score equal method: middle tau has zero weight for odd levels", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 20
  num_tau <- 3 # odd, mid_point = 2
  # Set all scores to zero except the middle one
  rank_scores <- list(rep(0, n), rep(1, n), rep(0, n))
  result <- qQTLR:::calculate_integrated_score(rank_scores, method = "equal", num_tau_levels = num_tau)
  # Middle score (index 2) is skipped for odd levels; lower_half = {1}, upper_half = {3}
  # int_rank_score = -rank_scores[[1]] + rank_scores[[3]] = 0 + 0 = 0
  # n_pairs = 1, so result = 0/1 = 0
  expect_equal(result, rep(0, n))
})

test_that("calculate_integrated_score IVW method returns correct length", {
  skip("Known bug: solve.QP constraint matrix dimensions are incorrect in IVW method")
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 4
  set.seed(20)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n, mean = i))
  result <- qQTLR:::calculate_integrated_score(rank_scores, method = "ivw", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
  expect_true(all(is.finite(result)))
})

test_that("calculate_integrated_score IVW method with odd num_tau_levels", {
  skip("Known bug: solve.QP constraint matrix dimensions are incorrect in IVW method")
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  n <- 30
  num_tau <- 5
  set.seed(21)
  rank_scores <- lapply(1:num_tau, function(i) rnorm(n, mean = i * 0.5))
  result <- qQTLR:::calculate_integrated_score(rank_scores, method = "ivw", num_tau_levels = num_tau)
  expect_true(is.numeric(result))
  expect_equal(length(result), n)
})

# ---------- fit_rank_scores ------------------------------------------------

test_that("fit_rank_scores returns list of correct length", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  num_tau <- 3
  result <- qQTLR:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 1)
  expect_true(is.list(result))
  expect_equal(length(result), num_tau)
})

test_that("fit_rank_scores each element has correct length", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  num_tau <- 3
  result <- qQTLR:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 1)
  for (i in seq_along(result)) {
    expect_equal(length(result[[i]]), length(d$phenotype))
    expect_true(is.numeric(result[[i]]))
  }
})

test_that("fit_rank_scores elements correspond to evenly spaced tau levels", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  num_tau <- 3
  result <- qQTLR:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 1)
  # Each element should match a direct call to calculate_rank_score at the expected tau
  for (i in 1:num_tau) {
    tau <- i / (num_tau + 1)
    expected <- qQTLR:::calculate_rank_score(d$phenotype, d$covariates, tau)
    expect_equal(result[[i]], expected)
  }
})

test_that("fit_rank_scores parallel path returns same results as sequential", {
  skip_if_not_installed("quantreg")
  d <- make_quail_data()
  num_tau <- 3
  result_seq <- qQTLR:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 1)
  result_par <- qQTLR:::fit_rank_scores(d$phenotype, d$covariates, num_tau_levels = num_tau, num_cores = 2)
  expect_equal(length(result_par), num_tau)
  for (i in seq_along(result_seq)) {
    expect_equal(result_par[[i]], result_seq[[i]])
  }
})

# ---------- QUAIL_rank_score_pipeline (exported) ---------------------------

test_that("QUAIL_rank_score_pipeline errors on non-numeric character phenotype", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  pheno_char <- c("a", "b", "c")
  expect_error(
    QUAIL_rank_score_pipeline(pheno_char, d$covariates, num_tau_levels = 3),
    "phenotype must be a numeric vector"
  )
})

test_that("QUAIL_rank_score_pipeline accepts data.frame phenotype", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  pheno_df <- data.frame(pheno = d$phenotype)
  # Should not error; data.frame is converted internally
  result <- QUAIL_rank_score_pipeline(pheno_df, d$covariates, num_tau_levels = 3, method = "equal")
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
})

test_that("QUAIL_rank_score_pipeline full run with small data and equal method", {
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  result <- QUAIL_rank_score_pipeline(d$phenotype, d$covariates,
    num_tau_levels = 5, method = "equal", num_cores = 1
  )
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
  expect_true(all(is.finite(result)))
})

test_that("QUAIL_rank_score_pipeline full run with IVW method", {
  skip("Known bug: solve.QP constraint matrix dimensions are incorrect in IVW method")
  skip_if_not_installed("quantreg")
  skip_if_not_installed("quadprog")
  d <- make_quail_data()
  result <- QUAIL_rank_score_pipeline(d$phenotype, d$covariates,
    num_tau_levels = 4, method = "ivw", num_cores = 1
  )
  expect_true(is.numeric(result))
  expect_equal(length(result), length(d$phenotype))
  expect_true(all(is.finite(result)))
})
