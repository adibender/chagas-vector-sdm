#' Add model prediction to data frame
#'
#' @param mod Model of choice.
#' @param newdata The data set to which prediction column will be added.
#' @param ... Further arguments passed to methods.
#' @importFrom dplyr mutate
#' @importFrom stats predict
#'
#' @export
add_prediction <- function(mod, newdata, ...) {
  UseMethod("add_prediction", mod)
}

#' @rdname add_prediction
#' @inherit add_prediction
#' @importFrom dplyr mutate
#' @export
add_prediction.gam <- function(mod, newdata, ...) {

  newdata %>%
    mutate(
      prediction = predict(mod, newdata, type = "response", ...))

}

#' @rdname add_prediction
#' @inherit add_prediction
#' @export
add_prediction.train <- function(mod, newdata, ...) {

   newdata %>%
    mutate(prediction = as.numeric(predict(mod, newdata, type = "prob")[, 2]))

}
