#' Binomial deviance evaluation metric for xgboost
#'
#' Can be specified as `eval_metric` in \code{\link[xgboost]{xgboost}}.
#' Calculates weighted binomial deviance
#'
#' @importFrom xgboost getinfo
#'
#' @export
xgb_binom_deviance <- function(preds, dtrain) {

  truth  <- getinfo(dtrain, "label")
  weight <- getinfo(dtrain, "weight")

  if(is.null(weight)) {
    weight <- rep(1, length(truth))
  }

  err <- binom_deviance(
    truth      = truth,
    prediction = preds,
    weights    = weight)

  list(
    metric = "binomial_deviance",
    value  = err)

}


#' @rdname eval_prev_gam
#' @inherit eval_prev_gam
#' @inherit as_xgb_matrix
#' @export
eval_prev_xgb <- function(mod, data, candidates,
  label = NULL, weight_col = NULL) {

  data_list <- data %>%
    select(id, one_of(c(candidates, label, weight_col)), no_tested, wght) %>%
    split_point_pol() %>%
    purrr::keep(~nrow(.x) > 0)
  pred_list <- data_list %>%
    map(~predict(
      mod,
      as_xgb_matrix(.x, mod$feature_names)))

  point_pol_eval(data_list, pred_list)

}

#' Transform data frame to xgb dense matrix
#'
#' Performs some preprocessing and creates a dense matrix for xgboost using
#' \code{\link[xgboost]{xgb.DMatrix}}.
#'
#' @param data A data frame
#' @param features A character vector with column names that will be kept in the
#' dense matrix.
#' @param label The name of the column containing the outcome
#' @param weights A (optional) numeric vector with sampling weights. Will be
#' added to the xbg dense matrix using \code{setinfo(..., "weight", weights)}
#' @param scale Logical. Should parameters be scaled before they are passed
#' to \code{\link[xgboost]{xgb.DMatrix}}
#'
#' @importFrom xgboost xgb.DMatrix
#' @importFrom xgboost setinfo
#'
#' @export
as_xgb_matrix <- function(
  data,
  features,
  label      = NULL,
  weight_col = NULL) {

  data <- na.omit(data[, c(features, label, weight_col)])


  form <- paste0("~", paste0(features, collapse = "+")) %>% as.formula()
  mm   <- model.matrix(form, data = data)[, -1]
  xgb_data <- xgb.DMatrix(data = mm)

  if (!is.null(label)) {
    setinfo(xgb_data, "label", as.numeric(data[[label]]))
  }

  if (!is.null(weight_col)) {
    setinfo(xgb_data, "weight", data[[weight_col]])
  }

  xgb_data

}


#' Run XGB cross validation given parameters
#'
#'
#' @param data Data frame containing all necessary variables.
#' @param xgb_params A list containing xgboost arguments. Will be used in
#' \code{do.call(xgboost, xgb_params)}.
#' @param features A character vector of covariates in \code{mod_data} to include
#' in model.
#' @param label The name of the outcome column.
#' @param deviance_base a list of baseline deviances in each fold to which
#' the deviance obtained by the xgb will be compared.
#' @param deviance_base_eval The same as above, but calculated on evaluation fold
#' instead of training folds.
#'
#' @importFrom xgboost xgboost
#' @importFrom purrr map_dfr
#' @importFrom dplyr filter
#' @importFrom tibble tibble
#'
#' @export
cv_iprev_xgb <- function(
  data,
  xgb_params,
  features,
  label,
  deviance_base,
  deviance_base_eval) {

  # run CV
  .folds             <- unique(data$fold) %>% sort()
  deviance_base      <- deviance_base[.folds]
  deviance_base_eval <- deviance_base_eval[.folds]
  map_dfr(
    .folds,
    function(.fold) {

      cat("Evaluation fold: ", .fold, "\n\n")

      # preparation
      train_df      <- dplyr::filter(data, fold != .fold)
      eval_df       <- dplyr::filter(data, fold  == .fold)
      xgb_params$data   <- as_xgb_matrix(train_df, features, label)
      xgb_params$weight <- train_df$wght2

      # fit model
      mod <- do.call(xgboost, xgb_params)

      # calculates deviance
      deviance_train <- eval_prev_xgb(mod, train_df, features, label, "wght2")
      deviance_eval <- eval_prev_xgb(mod, eval_df, features, label, "wght2")

      ## output
      xgb_params$data <- xgb_params$weight <- NULL
      tibble::tibble(
        eval_fold       = .fold,
        mod             = list(mod),
        xgb_params      = list(xgb_params),
        dev_train_point = deviance_train[1],
        dev_train_pol   = deviance_train[2],
        dev_eval_point  = deviance_eval[1],
        dev_eval_pol    = deviance_eval[2],
        dev_train_ratio_point = dev_train_point / deviance_base[[.fold]][1],
        dev_train_ratio_pol   = dev_train_pol / deviance_base[[.fold]][2],
        dev_eval_ratio_point  = dev_eval_point / deviance_base_eval[[.fold]][1],
        dev_eval_ratio_pol    = dev_eval_pol / deviance_base_eval[[.fold]][2])

    })

}
