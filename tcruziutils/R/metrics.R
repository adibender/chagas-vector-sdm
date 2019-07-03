#' Log-score
#'
#' @importFrom MLmetrics LogLoss
#' @rdname metrics
#' @export
log_loss <- function(truth, prediction) {
  LogLoss(prediction, truth)
}

#' Area under the curve
#' @importFrom MLmetrics AUC
#' @rdname metrics
#' @export
auc <- function(truth, prediction) {

  AUC(prediction, truth)

}

#' Brier score
#'
#'@rdname metrics
#' @export
brier <- function(truth, prediction) {
  mean( (truth - prediction) ** 2 )
}


#' Spiegelhalters z-Score
#'
#' Measure of model calibration (the lower the better)
#' @rdname metrics
#' @export
z_score <- function(truth, prediction) {

  nom <- sum( (truth - prediction) * (1 - 2 * prediction))
  denom <- sqrt(sum(
    ( (1 - 2 * prediction) ** 2) *
    (prediction * (1 - prediction)) ))

  abs(nom / denom)

}

#' @rdname metrics
#' @export
binom_deviance <- function(
  truth,
  prediction,
  weights = rep(1, length(truth))) {

  calc_deviance(truth, prediction, weights)

}

#' Caret train weighted deviance
#' @rdname metrics
#' @export
caret_deviance <- function(data, lev = NULL, model = NULL) {
  cls <- levels(data$obs) #find class names
  if (is.null(data$weights)){
    data$weights <- rep(1, nrow(data))
  }
  loss <- -binom_deviance(
    as.numeric(data$obs) - 1,
    data[, cls[2]],
    weights = data$weights)
  names(loss) <- c('BinomDeviance')
  loss
}


#' Caret train z-score function
#' @rdname metrics
#' @export
caret_z_score <- function(data, lev = NULL, model = NULL) {
  cls <- levels(data$obs) #find class names
  loss <- - z_score(
    as.numeric(data$obs) - 1,
    data[, cls[2]])
  names(loss) <- c('ZScore')
  loss
}

# from dismo::calc.deviance
calc_deviance <- function(obs, pred, weights = rep(1, length(obs)),
  family = "binomial", calc.mean = TRUE) {
    if (length(obs) != length(pred)) {
        stop("observations and predictions must be of equal length")
    }
    y_i <- obs
    u_i <- pred
    family = tolower(family)
    if (family == "binomial" | family == "bernoulli") {
        deviance.contribs <- (y_i * log(u_i)) + ((1 - y_i) *
            log(1 - u_i))
        deviance <- -2 * sum(deviance.contribs * weights)
    }
    else if (family == "poisson") {
        deviance.contribs <- ifelse(y_i == 0, 0, (y_i * log(y_i/u_i))) -
            (y_i - u_i)
        deviance <- 2 * sum(deviance.contribs * weights)
    }
    else if (family == "laplace") {
        deviance <- sum(abs(y_i - u_i))
    }
    else if (family == "gaussian") {
        deviance <- sum((y_i - u_i) * (y_i - u_i))
    }
    else {
        stop("unknown family, should be one of: \"binomial\", \"bernoulli\", \"poisson\", \"laplace\", \"gaussian\"")
    }
    if (calc.mean)
        deviance <- deviance/length(obs)
    return(deviance)
}


#' Weighted binomial deviance
#'
#' Calculates (weighted) binomial deviance for \code{mgcv::gam} models.
#'
#' @param mod The estimated GAM
#' @param test Test data (a data frame).
#'
#' @import dplyr
#' @importFrom stats predict
#' @importFrom mgcv predict.gam
#'
#' @export
eval_prev_gam <- function(mod, data, use_weight = TRUE) {

  data_list <- data %>% split_point_pol()

  data_list <- data_list %>% purrr::keep(~nrow(.x) > 0)

  pred_list <- data_list %>%
    map(~predict(mod, .x, type = "response", discrete = FALSE))

  point_pol_eval(data_list, pred_list, use_weight = use_weight)

}

split_point_pol <- function(data) {

   data <- data %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(n = dplyr::n())

  data_point   <- subset(data, n == 1) %>% as.data.frame()
  data_polygon <- subset(data, n > 1) %>% as.data.frame()

  data_list <-  list(point = data_point, polygon = data_polygon)

}


point_pol_eval <- function(data_list, pred_list, use_weight = TRUE) {

    map2_dbl(data_list, pred_list, ~{
    ind_valid <- !is.na(.y)
    if(use_weight) {
      binom_deviance(
        truth      = .x$prevalence[ind_valid],
        prediction = .y[ind_valid],
        weights    = .x$no_tested[ind_valid] * .x$wght[ind_valid])
    } else {
      binom_deviance(
        truth      = .x$prevalence[ind_valid],
        prediction = .y[ind_valid])
    }
  })

}
