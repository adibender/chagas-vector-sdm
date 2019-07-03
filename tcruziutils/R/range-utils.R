#' Add an element containing a variogram to a "fold" object
#'
#' @importFrom automap autofitVariogram
#' @importFrom stats glm binomial residuals
#' @export
get_variogram <- function(
  data,
  form = presence ~ elevation + rainfall + reflectance_evi + reflectance_tcb +
      reflectance_tcw + lst_day + lst_night + lst_diff) {

  # fit basic odel to presence data
  mod <- glm(form, data@data, family = binomial())

  # add residuals to training data
  resid_ind <- if (is.null(mod$na.action)) {
    1:nrow(data)
  } else {
    setdiff(1:nrow(data), mod$na.action)
  }
  data$resid <- NA
  data$resid[resid_ind] <- residuals(mod)

  # fit variogram to residuals
  automap::autofitVariogram(
    resid ~ 1,
    input_data = data[resid_ind, ],
    model = "Mat")

}
