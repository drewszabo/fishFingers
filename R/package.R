# R/package.R
#' @import data.table
#' @importFrom utils read.csv unzip
#' @importFrom rcdk parse.smiles get.fingerprint
#' @importFrom fingerprint fp.to.matrix
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols mutate
#' @importFrom stringr str_replace
#' @importFrom stats predict rbinom
#' @importFrom xgboost xgb.load xgb.DMatrix
#' @importFrom webchem is.smiles
#' @importFrom jsonlite fromJSON
#' @importFrom caret preProcess predict
#' @importFrom purrr map_dfr
NULL
