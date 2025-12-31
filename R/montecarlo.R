# R/montecarlo.R
#' Generate a single Monte Carlo binary matrix
#'
#' @param prob_matrix Matrix of probabilities (0–1)
#' @return Binary matrix with same dimensions and dimnames
generate_mc_matrix <- function(prob_matrix) {

  prob_matrix <- as.matrix(prob_matrix)

  nr <- nrow(prob_matrix)
  nc <- ncol(prob_matrix)

  mc_bin <- matrix(
    rbinom(nr * nc, size = 1, prob = prob_matrix),
    nrow = nr,
    ncol = nc,
    dimnames = dimnames(prob_matrix)
  )

  mc_bin
}


#' Generate deterministic binary fingerprints via Monte Carlo averaging
#'
#' @param prob_matrix Matrix of probabilities (0–1)
#' @param N Number of Monte Carlo samples
#' @param threshold Probability threshold for final binary call
#'
#' @return Deterministic binary matrix
generate_best_mc_binary <- function(prob_matrix,
                                    N = 10000,
                                    threshold = 0.5) {

  prob_matrix <- as.matrix(prob_matrix)

  nr <- nrow(prob_matrix)
  nc <- ncol(prob_matrix)

  avg_matrix <- matrix(0, nrow = nr, ncol = nc)

  for (i in seq_len(N)) {
    avg_matrix <- avg_matrix + generate_mc_matrix(prob_matrix)
  }

  avg_matrix <- avg_matrix / N

  final_binary <- (avg_matrix >= threshold) * 1

  dimnames(final_binary) <- dimnames(prob_matrix)

  final_binary
}


#' Monte Carlo binarisation for fishFingers output
#'
#' Automatically detects probabilistic fingerprint columns
#' by matching "Un" in column names.
#'
#' @param x Data frame returned by fishFingers()
#' @param N Number of Monte Carlo samples
#' @param threshold Binary threshold
#'
#' @return Data frame with Monte Carlo binary fingerprints
fishFingers_mc <- function(x,
                           N = 10000,
                           threshold = 0.5) {

  stopifnot(is.data.frame(x))

  # Identify probabilistic fingerprint columns
  prob_cols <- grep("Un", colnames(x), value = TRUE)

  if (length(prob_cols) == 0) {
    stop("No probability columns detected (no 'Un' in column names).")
  }

  prob_matrix <- as.matrix(
    as.data.frame(x[, ..prob_cols, drop = FALSE])
  )

  bin_matrix <- generate_best_mc_binary(
    prob_matrix = prob_matrix,
    N = N,
    threshold = threshold
  )

  bin_matrix
}
