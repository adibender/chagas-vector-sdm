#' Add weights and responses for poisson point process models
#'
#' @importFrom spatstat ppp quadscheme ppm
#' @importFrom dplyr mutate filter
#' @param ppm_data Presence only data padded with background points.
#' See \code{add_quads}.
#' @param spat_form RHS formula that indicates names of variables that
#' store longitude and latitude information (in that order).docu
#' @param presence_var The variable in \code{ppm_data} that stores information
#' on presence/absence.
#' @export
add_ppm_data <- function(
  ppm_data,
  presence_var = "presence",
  spat_form    = ~longitude + latitude) {

  lon_var <- all.vars(spat_form)[1]
  lat_var <- all.vars(spat_form)[2]

  win <- as_owin(ppm_data)
  presence_df <- filter(ppm_data, .data[[presence_var]] == 1)
  ppp <- ppp(presence_df[[lon_var]], presence_df[[lat_var]], window = win)
  absence_df <- filter(ppm_data, .data[[presence_var]] == 0)
  ppp_quad <- ppp(absence_df[[lon_var]], absence_df[[lat_var]], window = win)
  qs <- quadscheme(ppp, dummy = ppp_quad)

  wght <- ppm(qs, trend = ~1, preponly = TRUE)$Q$w
  ppm_data <- ppm_data %>%
    arrange(desc(.data[[presence_var]])) %>% # arrange b/c that how data is stored in Q$w
    mutate(wght = wght) %>%
    mutate(yw   = .data[[presence_var]] / .data$wght)

}

eval_mod <- function(
  mod,
  newdata,
  ...,
  truth = "presence",
  measures = c("auc", "brier", "logloss")) {

  prediction <- predict(mod, newdata, ...)
  truth <- newdata[[truth]]
  map(measures,
    ~do.call(.x, args = list(truth = truth, prediction = prediction)))
}

#' Transform object of class extent to object of class owin
#'
#'
#' @importFrom raster extent
#' @importFrom spatstat owin
#' @keywords internal
as_owin <- function(data, spat_form = ~longitude + latitude, poly = NULL) {

  extent <- data %>% as_spatial(spat_form) %>% extent()

  owin(
    xrange = c(extent@xmin, extent@xmax),
    yrange = c(extent@ymin, extent@ymax),
    poly = poly)

}

#' Wrapper function that performs PPM and creates maps
#'
#' @inheritParams get_pmap_ci
#' @inheritParams make_gam_formula
#' @param polygon A \code{SpatialPolygons*} object.
ppm_fun <- function(
  fold,
  candidates,
  env_grids,
  polygon,
  queries,
  gamma,
  ntile  = 200,
  method = "grid") {

  # basic quantities for species
  hull        <- fold$hull
  fold_extent <- raster::extent(hull)
  blocks      <- fold$blocks$blocks

  # data frame used for prediction
  pred_df <- env_grids %>%
    raster::crop(fold_extent) %>%
    raster::mask(hull) %>%
    grid_to_df()

  # clean data set (remove duplicates with respect to longigute, latitude,
    # merge_year, only keep variables in temp_form (as well as the above),
    # remove missing values
  temp_form <- paste0("~",
    paste0(c("id", "fold", candidates), collapse = "+")) %>%
    as.formula()
  train_df <- fold$train[fold$train$presence == 1, ] %>%
    remove_duplicates(form = temp_form)

  # add quadrature points for PPM
  # see ?spatstat::quadscheme for details
  ppm_data <- add_quads(
    data     = train_df,
    hull     = hull,
    queries  = queries,
    ntile    = ntile,
    year_var = "merge_year",
    method   = method)

  # add transformed outcome weights and offset
  ppm_data <- ppm_data %>% add_ppm_data()
  # add any fold, needed in make_gam_formula but not used as no
  # cross-validation will be performed for ppm models for weights
  ppm_data$fold <- 1
  # exclude locations in test data
  test_blocks <- blocks[blocks$folds == 5, ]
  ind <- sp::over(as_spatial(ppm_data), test_blocks)
  ppm_data <- list(
    ppm_train_df = ppm_data[!complete.cases(ind), ],
    ppm_test_df = ppm_data[complete.cases(ind), ])

  # visualization of training and test data
  p_bg_ppm <- tm_shape(fold$hull) + tm_borders(alpha = .1) +
    tm_shape(as_spatial(ppm_data[[1]][ppm_data[[1]]$presence == 0, ])) +
      tm_dots(col = "firebrick4", alpha = 0.5, size = .01) +
    tm_shape(as_spatial(fold$train[fold$train$presence == 1, ])) +
      tm_dots() +
    tm_shape(as_spatial(fold$test[fold$test$presence == 1, ])) +
      tm_dots() +
    tm_shape(test_blocks) +
      tm_borders(alpha = 1) +
    p_endemic_zone

  # fit the model
  ppm_form <- make_gam_formula(
    data = ppm_data$ppm_train_df[ppm_data$ppm_train_df$presence == 1, ],
    candidates = candidates,
    lhs = "yw~",
    type       = "smooth") %>%
  add_gp()
  ppm <- bam(
    ppm_form,
    data     = ppm_data$ppm_train_df,
    weights  = wght,
    family   = quasi(link = "log", variance = "mu"),
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1L)

  # pred_terms <- predict(ppm, pred_df, type = "terms", block.size = 5000,
  #   discrete = FALSE)
  pred_response <- predict(ppm, pred_df, type = "response",
    block.size = 1e4, discrete = FALSE)
  pred_df$prediction <- as.numeric(pred_response)
  pred_grid <- df_to_grid(pred_df, env_grids[[1]], "prediction") %>%
    raster::crop(fold_extent)
  max_pred <- max(raster::values(pred_grid), na.rm = TRUE)
  # raster::values(pred_grid) <- raster::values(pred_grid) * 5
  # brks <- seq(0, 1, by = .2)
  brks <- c(0, exp(seq(0, log(max_pred), length.out = 4)))

  p_pred <-
    tm_shape(raster::crop(countries, fold_extent), is.master = TRUE) +
      tm_borders(alpha = .5) +
    tm_shape(fold$hull) + tm_borders(alpha = .5) +
    tm_shape(pred_grid) +
      tm_raster(alpha = .8,
        breaks = brks,
        title = "prediction",
        palette = viridis::magma(1e3),
        style = "cont") +
    tm_layout(
      title             = fold$species,
      title.size        = 1.5,
      title.fontface    = "italic",
      legend.position   = c("left", "bottom"),
      legend.text.size  = 1.2,
      legend.hist.size  = 1.2,
      legend.title.size = 1.5,
      bg.color          = "whitesmoke")

  list(species = fold$species, model = ppm,
    p1 = p_pred, p2 = p_pred2, p_gp = p_gp, p_bg = p_bg_ppm)

}
