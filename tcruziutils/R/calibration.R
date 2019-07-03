#' Add calibration data to df
#'
#' @importFrom dplyr mutate select group_by n
#' @param data A data set containing the true (binary) outcome and
#' a column with predictions.
#' @param truth Name of the column containing true outcome
#' @param ncut The number of intervals in which predictions will be disretized.
#' For each interval, the percentage of events (outcome = 1) is calculated.
#' @param metrics A named vector of functions (metrics) that will be applied
#' to outcomes and predictions.
#' @export
calc_calibration <- function(
  data,
  truth = "presence",
  ncut      = 50L,
  metrics   = c(
    auc      = auc,
    brier    = brier,
    z_score  = z_score,
    deviance = binom_deviance)) {

  ndf <-  data %>%
    select(prediction, one_of(truth)) %>%
    na.omit() %>%
    mutate(bin.pred = cut(.data$prediction, ncut)) %>%
    group_by(bin.pred) %>%
    summarize(
      bin.y = mean(.data[[truth]]),
      bin.yhat = mean(prediction),
      n = n()) %>%
    ungroup() %>%
    mutate(cal_curve = fitted(gam(bin.y ~ s(bin.yhat, bs = "ts"))))

  scores <- scores(data)
  attr(ndf, "scores") <- scores
  print(scores)

  return(ndf)

}

#' @rdname calc_calibration
#' @inherit calc_calibration
#' @export
calc_calibration_sp <- function(
  data,
  truth = "presence",
  ncut      = 50L,
  metrics   = c(
    auc      = auc,
    brier    = brier,
    z_score  = z_score,
    deviance = binom_deviance)) {

  ndf <-  data %>%
    select(layer, prediction, one_of(truth)) %>%
    na.omit() %>%
    group_by(layer) %>%
    summarize(
      bin.y = mean(.data[[truth]]),
      bin.yhat = mean(prediction),
      n = n()) %>%
    ungroup() %>%
    mutate(cal_curve = fitted(gam(bin.y ~ s(bin.yhat, bs = "ts"))))

  scores <- scores(data)
  attr(ndf, "scores") <- scores
  print(scores)

  return(ndf)

}

#' Plot calibration curve
#'
#' @import ggplot2
#' @param calibration A data frame containing calibration info.
#' See \code{\link{calc_calibration}}
#' @export
gg_calibration <- function(calibration) {

  gg_cal <- ggplot(calibration, aes(x = bin.yhat, y = bin.y)) +
    geom_point(aes(size = n)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    geom_abline(intercept = 0, slope = 1, col = 1, linetype = 2) +
    # geom_line(aes(y = cal_curve), linetype = 2) +
    geom_smooth(se = TRUE, method="gam", formula = y~s(x, bs = "ts")) +
    geom_smooth(method = "lm", se = TRUE, col = "brown")

  gg_cal

}

#' Calculate various classification metrics
#'
#' @rdname metrics
#' @importFrom dplyr select one_of mutate_if rename
#' @importFrom purrr map_dbl
#' @export
scores <- function(
  data,
  pred_col  = "prediction",
  truth = "presence",
  metrics   = c(
    auc      = auc,
    brier    = brier,
    z_score  = z_score,
    deviance = binom_deviance)) {

  data <- data %>%
    select(one_of(c(truth, pred_col))) %>%
    mutate_if(is.factor, ~as.numeric(.)-1) %>%
    rename(truth = truth, prediction = pred_col) %>%
    na.omit()
  map_dbl(metrics, ~ exec(.x, !!!as.list(data)))

}
