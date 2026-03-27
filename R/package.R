#' R/package.R
#' @import data.table
#' @importFrom utils read.csv unzip
#' @importFrom rcdk parse.smiles get.fingerprint
#' @importFrom fingerprint fp.to.matrix
#' @importFrom tibble as_tibble tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_cols mutate
#' @importFrom stringr str_replace
#' @importFrom stats predict rbinom
#' @importFrom xgboost xgb.load xgb.DMatrix
#' @importFrom webchem is.smiles
#' @importFrom jsonlite fromJSON
#' @importFrom caret preProcess
#' @importFrom purrr map_dfr map_df
#' @importFrom RSirius SiriusSDK
NULL
