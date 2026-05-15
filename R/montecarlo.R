# R/montecarlo.R

#' Generate deterministic binary fingerprints via vectorized Binomial sampling
#'
#' @param prob_matrix Matrix of probabilities (0–1)
#' @param N Number of Monte Carlo samples
#' @param threshold Probability threshold for final binary call
#' @return Deterministic binary matrix
generate_best_mc_binary <- function(prob_matrix, N = 10000, threshold = 0.5) {
  
  prob_matrix <- as.matrix(prob_matrix)
  
  # 1. Validation (Done once, not N times)
  if (!is.numeric(prob_matrix)) {
    storage.mode(prob_matrix) <- "numeric"
  }
  if (any(is.na(prob_matrix))) stop("prob_matrix contains NA values.")
  if (any(prob_matrix < 0 | prob_matrix > 1)) stop("Values must be between 0 and 1.")

  # 2. Vectorized Math Replacement
  # Instead of looping N times, we use the property that the sum of N 
  # Bernoulli trials is a Binomial distribution with size = N.
  sum_matrix <- matrix(
    rbinom(length(prob_matrix), size = N, prob = prob_matrix),
    nrow = nrow(prob_matrix),
    ncol = ncol(prob_matrix)
  )

  # 3. Average and Threshold
  final_binary <- ((sum_matrix / N) >= threshold) * 1
  dimnames(final_binary) <- dimnames(prob_matrix)
  
  return(final_binary)
}

#' Monte Carlo binarisation for fishFingers output
#'
#' @param x Data frame returned by fishFingers()
#' @param N Number of Monte Carlo samples
#' @param threshold Binary threshold
#' @return Data frame with Monte Carlo binary fingerprints
fishFingers_mc <- function(x, N = 10000, threshold = 0.5) {
  stopifnot(is.data.frame(x))

  prob_cols <- grep("Un", colnames(x), value = TRUE)
  if (length(prob_cols) == 0) stop("No probability columns detected.")

  prob_matrix <- as.matrix(x[, prob_cols, drop = FALSE])

  # Calls the optimized function above
  bin_matrix <- generate_best_mc_binary(
    prob_matrix = prob_matrix,
    N = N,
    threshold = threshold
  )

  return(bin_matrix)
}