#' Create thematic map, possibly stratified by year
#'
#' A high level plotting function to gain quick overview over spatial
#' distribution of environmental variables, possibly stratified by year.
#'
#' @import tmap
#' @importFrom purrr imap
#' @importFrom furrr future_imap
#' @importFrom mastergrids read_grids
#'
#' @param grids A list of gridded data sets (SpatialGridDataFrame). Each
#' list entry corresponds to observations from one year.
#' @inheritParams tmap::tm_raster
#' @inheritParams tmap::tm_layout
#' @inheritParams tmap::tmap_arrange
#'
#' @export
tmap_grid_by_year <- function(
  grids,
  palette         = NULL,
  breaks          = NULL,
  legend.position = c("left", "bottom"),
  ncol            = NA,
  style           = "cont",
  midpoint        = NULL,
  alpha_raster    = 0.8,
  outer.margins   = 0,
  base_plot = tcruziutils::p_endemic_zone,
  ...) {

  if (is.character(grids)) {
    grid_names <- grids
    grids <- read_grids(grids,  ...)
    names(grids) <- grid_names
  }
  tmap_list <- future_imap(grids,
    ~ p_endemic_zone +
      tm_shape(.x) +
        tm_raster(col = names(.x@data)[1], alpha = alpha_raster, palette =
         palette, breaks = breaks, title = "", style = style,
         midpoint = midpoint) +
      tm_layout(
        title = paste0("year = ", str_extract_year(.y)),
        legend.position = legend.position)) %>%
  c(outer.margins = outer.margins, ncol = ncol)
do.call(tmap_arrange, tmap_list)

}


#' Plot spatial distribution of data with different geometries
#'
#' Plots polygon data on administrative districts levels 1 and 2 together with
#' point data.
#'
#' @import tmap
#' @importFrom dplyr filter
#' @importFrom sp coordinates
#'
#' @param data A data set with variables Longitude and Latitude storing
#' spatial coordinates. Can be specified via \code{spat_form} argument.
#' @param shp_admin1 A \code{SpatialPolygonsDataFrame} with administrative
#' districts polygons level 1.
#' @param shp_admin2 As \code{shp_admin1} but for administrative districts
#' level 2.
#' @param keep The covariate of interest that will be used as intensity
#' or fill color.
#' @param spat_form Specifies longitude and latitude variables in
#' \code{data}. See \code{\link[sp]{coordinates}}.
#' @param base_plot A base plot (here endemic zone in south america) to which
#' polygon data and  point data will be added as additional layers.
#'
#' @export
spat_plot <- function(
  data,
  shp_admin1,
  shp_admin2,
  keep = "prevalence",
  spat_form = ~ longitude + latitude,
  base_plot = p_endemic_zone) {

  if(keep == "prevalence") {
    palette_admin1 = "Greens"
    palette_admin2 = "Blues"
    palette_dots   = "cividis"
  } else {
    palette_admin1 = "Set3"
    palette_admin2 = "Set2"
    palette_dots   = "Set1"
  }

  lon_var <- all.vars(spat_form)[1]
  data_lonlat <- data %>% filter(!is.na(.data[[lon_var]]))
  coordinates(data_lonlat) <- spat_form

  if ("admin1_1" %in% colnames(data)) {
    data_admin1 <- merge_admin(data, shp_admin1, keep = keep,
      by.y = "admin1_1")
  }
  if ("admin2_1" %in% colnames(data)) {
      data_admin2 <- merge_admin(data, shp_admin2, keep = keep)
  }

  if (exists("data_admin1")) {
    base_plot <- base_plot + tm_shape(data_admin1) +
      tm_polygons(col = keep, title = "Admin1", palette = palette_admin1)
  }
  if (exists("data_admin2")) {
    base_plot <- base_plot + tm_shape(data_admin2) +
      tm_polygons(col = keep, title = "Admin 2", palette = palette_admin2, lwd = 0.2)
  }

  base_plot +
    tm_shape(data_lonlat) + tm_dots(size = .05, alpha = 0.6, col = keep,
      title = "Coordinates", palette = palette_dots, shape = 19)
}

#' @inherit spat_plot
#' @export
#' @keywords internal
spat_plot_prep <- function(
  data,
  shp_admin1,
  shp_admin2,
  keep = "prevalence",
  spat_form = ~ longitude + latitude) {

  lon_var <- all.vars(spat_form)[1]
  data_lonlat <- data %>% filter(!is.na(.data[[lon_var]]))
  coordinates(data_lonlat) <- spat_form

  if ("admin1_1" %in% colnames(data)) {
    data_admin1 <- merge_admin(data, shp_admin1, keep = keep,
      by.y = "admin1_1")
  }
  if ("admin2_1" %in% colnames(data)) {
      data_admin2 <- merge_admin(data, shp_admin2, keep = keep)
  }

  list(points = data_lonlat, admin1_df = data_admin1, admin2_df = data_admin2)

}


#' Thematic maps of species presence
#'
#'
#' @export
tmap_species <- function(p_endemic_zone, folds, spec) {

  fold <- folds[[spec]]
  p_endemic_zone +
    tm_shape(fold$hull) + tm_borders() +
    tm_shape(rbind(fold$train, fold$test)) +
    tm_dots(col = "presence", style = "cat", palette = c("steelblue", "black"),
      alpha = .9)

}

#' @rdname tmap_species
#' @export
tmap_presence <- function(p_endemic_zone, folds, spec) {

  fold <- folds[[spec]]
  p_endemic_zone +
    tm_shape(fold$hull) + tm_borders() +
    tm_shape(subset(fold$train, presence == 1)) + tm_dots() +
    tm_shape(subset(fold$test, presence == 1)) +
      tm_dots(col = "firebrick4")

}

#' @rdname tmap_species
#' @export
tm_spec_pred <- function(p_endemic_zone, folds, spec,
  path = "prediction-maps/") {

  pres_map <- tmap_species(p_endemic_zone, folds, spec)
  pred_grid <- raster::raster(paste0(pred_files,
    sub(" ", "-", spec), "-prediction.tif"))

  pres_map +
    tm_shape(pred_grid) +
      tm_raster(
        alpha = .7,
        style = "cont",
        breaks = seq(0, 1, by = .2),
        palette = viridis::magma(1e3)) +
    tm_layout(
      title             = spec,
      title.size        = 1.5,
      title.fontface    = "italic",
      legend.position   = c("left", "bottom"),
      legend.text.size  = 1.2,
      legend.hist.size  = 1.2,
      legend.title.size = 1.5,
      bg.color          = "whitesmoke")


}

#' @rdname tmap_species
#' @export
tm_spec_pred_cut <- function(p_endemic_zone, folds, spec, q = .95) {

  pres_map <- tmap_species(p_endemic_zone, folds, spec)
  pred_grid <- raster::raster(paste0(pred_files,
    sub(" ", "-", spec), "-prediction.tif"))
  raster::values(pred_grid) <- 1L*(raster::values(pred_grid) >= .95)

  pres_map +
    tm_shape(pred_grid) +
      tm_raster(
        alpha = .7,
        style = "cat",
        title = "presence (95% quantile)") +
    tm_layout(
      title             = spec,
      title.size        = 1.5,
      title.fontface    = "italic",
      legend.position   = c("left", "bottom"),
      legend.text.size  = 1.2,
      legend.hist.size  = 1.2,
      legend.title.size = 1.5,
      bg.color          = "whitesmoke")


}

#' @export
tmap_cv <- function(
  fold,
  countries,
  palette = RColorBrewer::brewer.pal(9, "Set1")[1:2]) {

  species <- fold$species

  new_hull <- fold$train %>%
    rgeos::gConvexHull() %>%
    rgeos::gBuffer(width = 10)
  block <- fold$blocks$blocks

  mask <- fold$hull

  folds_train <- unique(fold$train$fold)
  folds_test <- unique(fold$test$fold)
  crp_countries <- raster::crop(countries, raster::extent(new_hull))
  p1 <- tm_shape(crp_countries, is.master = TRUE) +
    tm_borders(alpha = .5) +
    tm_shape(fold$hull) + tm_borders(alpha = .5) +
    tm_shape(mask) + tm_borders(alpha = 0.8) +
    tm_shape(fold$train[fold$train$presence == 0, ]) +
      tm_dots(
        col = "presence", style = "cat", palette = palette[2],
        alpha = 0.8, legend.show = FALSE) +
    tm_shape(fold$test[fold$test$presence == 0, ]) +
      tm_dots(
        col = "presence", style = "cat", palette = palette[2],
        alpha = 0.8, legend.show = FALSE) +
    tm_shape(fold$train[fold$train$presence == 1, ]) +
      tm_dots(
        col = "presence", style = "cat", palette = palette[1],
        alpha = 0.8, legend.show = FALSE) +
    tm_shape(fold$test[fold$test$presence == 1, ]) +
      tm_dots(
        col = "presence", style = "cat", palette = palette[1],
        alpha = 0.8, legend.show = FALSE) +
    tm_add_legend("symbol", col = palette[2:1], size = c(1, 1),
      labels = c("background", "presence"))
  p2 <- p1 +
    tm_shape(block[block$folds %in% folds_train, ]) +
      tm_borders() +
      tm_text(text = "folds", palette = "black", size = 0.7) +
    tm_shape(block[block$folds %in% folds_test, ]) +
      tm_polygons(palette = "grey70", alpha = 0.7) +
      tm_text(text = "folds", palette = "black", size = 0.7) +
    tm_layout(
      title          = species,
      title.fontface = "italic",
      title.size     = 1.5)

  p2

}

#' Map of best model selected via cross validation
#'
#'
#' @inheritParams get_pmap_ci
#' @export
tmap_sdm <- function(cv, fold, countries, width = 10) {

  # print(cv$species)
  tm_shape(countries, is.master = TRUE) +
    tm_borders(alpha = .5) +
  tm_shape(fold$hull) + tm_borders(alpha = .5) +
  tm_shape(cv$grid) +
    tm_raster(alpha = .7, breaks = seq(0, 1, by = .2),
      title = "prediction",
      # title = paste0("prediction (AUC: ", get_test_auc(cv), ")"),
      palette = viridis::magma(1e3), style = "cont") +
  tm_layout(
    title             = cv$species,
    title.size        = 1.5,
    title.fontface = "italic",
    legend.position   = c("left", "bottom"),
    legend.text.size  = 1.2,
    legend.hist.size  = 1.2,
    legend.title.size = 1.5,
    bg.color          = "whitesmoke")

}


#' Map with three panels
#'
#' Creates figure with three panels, lower and upper CI as well
#' as prediction (expectation).
#'
#' @param cv An object returned by \code{\link{cv_settings}}.
#' @param fold An object returned by \code{\link{get_sp_folds}}.
#' @param env_grids A \code{RasterBrick} object containing all
#' covariate layers used at estimation stage. Maps will reflect
#' model prediction for the area defined by \code{env_grids}.
#' @param countries A \code{SpatialPolygons*} object. Used as
#' background for the grid layer.
#' @param width Extra buffer around the extended hull used for
#' modelling for which countrie borders will also be shown.
#' Defaults to 10 degrees.
#' @importFrom mastergrids grid_to_df df_to_grid
#' @export
get_pmap_ci <- function(
  cv,
  fold,
  env_grids,
  countries,
  width = 10,
  discrete = FALSE) {

  new_hull <- fold$train %>%
    rgeos::gConvexHull() %>%
    rgeos::gBuffer(width = 10)
  crp_countries <- raster::crop(countries, raster::extent(new_hull))

  ndf <- env_grids %>%
    raster::crop(raster::extent(fold$hull)) %>%
    raster::mask(fold$hull) %>%
    mastergrids::grid_to_df()

  pred <- predict(get_best_mod(cv), newdata = ndf, type = "link", se.fit = TRUE,
    discrete = discrete, block.size = 1e5)
  # fit <- exp(pred$fit)/(1 + exp(pred$fit))
  ci_lower <- pred$fit - 2*pred$se
  ci_lower <- exp(ci_lower)/(1 + exp(ci_lower))
  ci_upper <- pred$fit + 2*pred$se
  ci_upper <- exp(ci_upper)/(1 + exp(ci_upper))

  # ndf$prediction = fit
  ndf$ci_lower <- ci_lower
  ndf$ci_upper <- ci_upper

  # grid_pred <- mastergrids::df_to_grid(ndf, env_grids[[1]], "prediction")
  grid_lower <- mastergrids::df_to_grid(ndf, env_grids[[1]], "ci_lower")
  grid_upper <- mastergrids::df_to_grid(ndf, env_grids[[1]], "ci_upper")

  map_pred <- tmap_sdm(cv, fold, crp_countries, width = 10)
  cv$grid <- grid_lower
  map_lower <- tmap_sdm(cv, fold, crp_countries, width = 10)
  cv$grid <- grid_upper
  map_upper <- tmap_sdm(cv, fold, crp_countries, width = 10)

  list(lower = map_lower, pred = map_pred, upper = map_upper)

}

#' Creates a "breaks" sequence for plotting
#'
#'
#' @export
#' @param vec A numeric vector.
get_breaks <- function(vec) {
  c(0, exp(seq(0, log(max(vec, na.rm = TRUE)), length.out = 5L)))
}


#' Bivariate map
#'
#' Bivariate map that simultaneously conveys prediction and uncertainty
#'
#' @param data A data frame containing colums "prediction" and "se".
#' @param grid A place-holder grid/raster.
#' @param fold Species information.
#' @param palette A bivariate colour palette.
#'
#' @import pals ggplot2
#' @importFrom classInt classIntervals findCols
#' @importFrom mastergrids df_to_grid grid_to_df
#' @importFrom sf st_bbox st_as_sf st_sfc
#' @importFrom stars st_as_stars geom_stars
#' @importFrom rgeos gBuffer
#' @importFrom tibble tibble
#' @importFrom raster crop projection
#' @export
tm_bivariate <- function(
  data,
  grid,
  fold,
  x_name = "Standard Error",
  y_name = "Prediction",
  palette    = stevens.bluered(n = 9),
  alpha      = .8,
  main.title = NULL,
  shp        = countries,
  size_frac  = 15,
  name       = "prediction") {

  hull <- fold$hull
  ehull <- fold$train %>%
    rgeos::gConvexHull() %>%
    rgeos::gBuffer(width = 10)

  n_brks <- sqrt(length(palette))
  brks_fit <- classInt::classIntervals(data$prediction, n = n_brks)
  brks_se <- classInt::classIntervals(data$se, n = n_brks)

  class_fit <- classInt::findCols(brks_fit)
  class_se <- classInt::findCols(brks_se)
  data$pred_c2 <- class_fit + n_brks * (class_se - 1)
  unique_val <- unique(data$pred_c2) %>% sort()
  pal_order <- map(seq_len(n_brks), ~.x + n_brks * (seq_len(n_brks) - 1)) %>%
    unlist()

  grid <- df_to_grid(data, grid, "pred_c2")
  names(grid) <- name

  shp <- shp %>% raster::crop(raster::extent(ehull))

  bbox_ehull <- st_bbox(shp)
  bbox_hull <- st_bbox(hull)
  star <- st_as_stars(
    grid,
    xlim = bbox_hull[c(1, 3)],
    ylim = bbox_hull[c(2, 4)])
  star[[name]] <- factor(star[[name]], levels = seq_along(palette))

  ## ggplot map object
  gg_map <- ggplot() +
    geom_sf(data = st_as_sf(shp), fill = NA, alpha = alpha) +
    coord_sf(
      default = TRUE,
      xlim    = bbox_ehull[c(1, 3)],
      ylim    = bbox_ehull[c(2, 4)]) +
    theme(panel.grid.major = element_line(colour = 'transparent'))
  ## add bivariate map raster layer
  gg_map <- gg_map +
    stars::geom_stars(
      data = star,
      mapping = aes_string(x = "x", y = "y", fill = name),
      alpha = alpha) +
    scale_fill_manual(values = palette[pal_order[unique_val]], guide = FALSE) +
    xlab("") + ylab("") + ggtitle(main.title) +
    theme(panel.grid.major = element_line(colour = 'transparent'))

  ## ggplot legend object
  legend <- tm_bivar_legend(
    palette         = palette,
    hull            = ehull,
    size_frac       = size_frac)
  gg_legend <- ggplot() +
    geom_sf(data = legend, mapping = aes(fill = fill)) +
    scale_fill_identity() +
    labs(
      x = bquote(.(x_name) ~ ""%->%""),
      y = bquote(.(y_name) ~ ""%->%"")) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = rel(1)),
    panel.grid.major = element_line(colour = 'transparent'))
  # print(gg_map)
  # vp <- grid::viewport(x = .2, y = .2, width = .2, height = .2)
  # grid::pushViewport(vp)
  # print(gg_legend, newpage = FALSE)
  # grid::popViewport()
  list(map = gg_map, legend = gg_legend)

}


#' @rdname tm_bivariate
#' @keywords internal
#' @importFrom sf st_sf st_polygon
#' @importFrom purrr map
#' @importFrom sp CRS
tm_bivar_legend <- function(
  palette,
  hull,
  size_frac = 15) {

  crs_hull <- CRS(projection(hull))
  bbox <- st_bbox(hull)

  lx <- bbox[3] - bbox[1]
  ly <- bbox[4] - bbox[2]
  step_length <- max(lx, ly) / 15
  b1 <- bbox[1]
  b2 <- bbox[2]

  first_rectangle <- c(
    b1 + 1 * step_length, b2 + 1 * step_length,
    b1 + 2 * step_length, b2 + 1 * step_length,
    b1 + 2 * step_length, b2 + 2 * step_length,
    b1 + 1 * step_length, b2 + 2 * step_length,
    b1 + 1 * step_length, b2 + 1 * step_length)
  # Copy it 8 times
  list_of_rectangles <- list(
    first_rectangle + c(0 * step_length, 0 * step_length),
    first_rectangle + c(1 * step_length, 0 * step_length),
    first_rectangle + c(2 * step_length, 0 * step_length),
    first_rectangle + c(0 * step_length, 1 * step_length),
    first_rectangle + c(1 * step_length, 1 * step_length),
    first_rectangle + c(2 * step_length, 1 * step_length),
    first_rectangle + c(0 * step_length, 2 * step_length),
    first_rectangle + c(1 * step_length, 2 * step_length),
    first_rectangle + c(2 * step_length, 2 * step_length)
  )

  # Transform it as a list of polygons
  geo_rectangles <- map(
    map(list_of_rectangles, ~ list(matrix(data = .x, ncol = 2, byrow = TRUE))),
    st_polygon)

    # Create the sf
    st_sf(
      tibble(fill = palette),
      geometry = geo_rectangles %>% st_sfc(),
      crs = crs_hull)
}


#' @rdname tm_bivariate
#' @param cv Cross validation output for species.
#' @param fold Fold information for species.
#' @param env_grids Environmental covariates
#' @param ... Further arguments passed to tm_bivariate
#' @importFrom mastergrids grid_to_df
#' @export
tm_bivar_cv <- function(cv, fold, env_grids, discrete = TRUE, ...) {

  mod <- get_best_mod(cv)
  hull <- fold$hull

  pred_df <- grid_to_df(env_grids, hull, hull)
  pred <- predict(get_best_mod(cv), pred_df, type = "response", se = TRUE,
    discrete = discrete,
    block.size = 1e5)
  pred_df$prediction <- pred$fit
  pred_df$se <- pred$se.fit

  tm_bivariate(pred_df, env_grids[[1]], fold, ...)

}

tm_bivar_draw <- function(tm_bivar,
  x = .15, y = .15, width = .25, height = .25) {

 cowplot::ggdraw() +
  cowplot::draw_plot(tm_bivar[[1]] +
    theme(axis.text = element_text(size = rel(1.2)))) +
  cowplot::draw_plot(tm_bivar[[2]], x, y, width, height)

}

#' Quick plot of predictions
#'
#' @export
tmap_pred <- function(pred_df, grid, shp, col = "prediction") {

  names(grid) <- col
  pred_grid <- mastergrids::df_to_grid(pred_df, grid, col)
  tm_shape(shp)  +
    tm_borders() +
  tm_shape(pred_grid) +
    tm_raster(alpha = .8, style = "cont", palette = viridis::magma(1e3))

}
