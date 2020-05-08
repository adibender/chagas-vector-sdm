#' Add background data and calculate spatial folds
#'
#' @importFrom blockCV spatialBlock
#' @importFrom raster crop extent raster area
#' @importFrom sp coordinates proj4string
#' @importFrom rgeos gBuffer gConvexHull
#'
#' @inheritParams blockCV::spatialBlock
#' @inheritParams rgeos::gBuffer
#' @param data A data set with spatial coordinates (but not a \code{Spatial*}
#' data set).
#' @param species The name of the species for which data will be generated.
#' @param species_var The name of the column in which names of species are
#' stored.
#' @param n_blocks The approximate number of blocks with which the spatial
#' region will be filled with. Each fold than consists of approximately
#' \code{n_blocks}/\code{k} blocks.
#' @param test_fold Indicates the block numbers that will be associated with
#' the test data set. Can be a vector. Defaults to the last fold, i.e.
#' blocks numbered \code{k}.
#' @param width The width by which convex hull around presence locations is
#' extended (in degrees).
#' @param ... Further arguments passed to \code{\link[blockCV]{spatialBlock}}
#' @export
get_sp_folds <- function(
  data,
  species,
  mask,
  rasterLayer   = NULL,
  k             = 5,
  n_blocks      = 200,
  test_fold     = k,
  species_var   = "species2",
  coord_form    = ~ longitude + latitude,
  width         = 5, # degrees
  extent        = NULL,
  proj          = proj_tcruzi,
  selection     = "random",
  showBlocks    = FALSE,
  progress      = FALSE,
  maskBySpecies = TRUE,
  calc_range    = TRUE,
  ... ) {

  data$presence <- 1L * (data[[species_var]] == species)
  if (!grepl("^Spatial.*", class(data)[1])) {
    coordinates(data) <- coord_form
  }
  data@data <- as.data.frame(data@data)
  if (!is.null(mask)) {
    proj <- proj4string(mask)
    extent <- extent(mask)
  } else {
    if (!is.null(rasterLayer)) {
      proj <- proj4string(rasterLayer)
      extent <- extent(rasterLayer)
    }
  }
  proj4string(data) <- proj
  hull   <- gConvexHull(data[data$presence == 1, ])
  proj4string(hull) <- ""
  hull <- rgeos::gBuffer(hull, width = width)
  proj4string(hull) <- proj
  hull   <- crop(hull, extent(mask))
  extent <- extent(hull)
  if (!is.null(rasterLayer)) {
    if (class(rasterLayer) == "SpatialGridDataFrame") {
      rasterLayer <- raster(rasterLayer)
    }
    rasterLayer <- crop(rasterLayer, extent)
  }
  if (!is.null(mask)) {
    mhull <- raster::intersect(hull, mask)
  }
  area <- area(mhull, byid = FALSE) %>% sum()
  # exclude absences outside extended hull
  data <- raster::intersect(data, hull)

  if (calc_range) {
    vario <- get_variogram(data)
    range <- vario$var_model$range[2] * 1e3
  } else {
    vario = NULL
    range = NULL
  }

  # approximates the length of one block based on the landmass area
  # of the spatial extent of the species and approximate number of blocks
  theRange <- sqrt(area / n_blocks)

  blocks <- spatialBlock(
    speciesData   = data,
    species       = "presence",
    rasterLayer   = rasterLayer,
    theRange      = theRange,
    k             = k,
    selection     = selection,
    biomod2Format = FALSE,
    showBlocks    = showBlocks,
    progress      = progress,
    maskBySpecies = maskBySpecies,
    ...)

  # add fold information to data set
  data$fold <- blocks$foldID
  data      <- data[!is.na(data$fold), ]

  list(
    species     = species,
    hull        = mhull,
    blocks      = blocks,
    extent      = extent,
    train       = data[!(data$fold %in% test_fold), ],
    test        = data[data$fold %in% test_fold, ],
    block_width = theRange,
    vario       = vario,
    vario_range = range)

}

#' Visualize spatial CV blocks
#'
#' Given an object returned by \code{get_sp_folds}, creates thematic maps
#' that visualize the spatial distribution of presence and absence points
#' as well as the location of blocks associated with different folds.
#'
#' @param blocks An object returned by \code{get_sp_folds}
#' @param mask An object of class \code{SpatialPointsPolygon} (or equivalent
#' object (e.g. \code{sf})) defining the spatial region in which observations
#' occurrred
#' @param palette_dots A vector of colors for absence and presence points
#' respectively.
#'
#' @importFrom tmap tm_shape tm_dots tm_borders tm_polygons tm_text tmap_arrange
#' @importFrom raster crop
#' @export
tmap_blocks <- function(
  blocks,
  mask         = NULL,
  palette_dots = c("steelblue", "black")) {

  if (is.null(mask)) {
    mask <- blocks$hull
  } else {
    mask <- crop(mask, blocks$extent)
  }
  # block polygon objects
  block       <- blocks$blocks$blocks
  folds_train <- unique(blocks$train$fold)
  folds_test  <- unique(blocks$test$fold)

  # make tmaps + combine
  p1 <- tm_shape(mask) + tm_borders(alpha = .8) +
    tm_shape(blocks$train[blocks$train$presence == 0, ]) +
      tm_dots(col = "presence", style = "cat", palette = palette_dots[1],
        alpha = .7, size = .075, legend.show = FALSE) +
    tm_shape(blocks$test[blocks$test$presence == 0, ]) +
      tm_dots(col = "presence", style = "cat",palette = palette_dots[1],
        alpha = .7, size = .075, legend.show = FALSE) +
    tm_shape(blocks$train[blocks$train$presence == 1, ]) +
      tm_dots(col = "presence", style = "cat", palette = palette_dots[2],
        alpha = .7, size = .075, legend.show = FALSE) +
    tm_shape(blocks$test[blocks$test$presence == 1, ]) +
      tm_dots(col = "presence", style = "cat",palette = palette_dots[2],
        alpha = .7, size = .075, legend.show = FALSE) +
    tm_add_legend("symbol", col = palette_dots, size = c(1, 1),
      labels = c("background", "presence")) +
    tm_layout(
      title          = blocks$species,
      title.size     = 1.5,
      title.fontface = "italic")
    p2 <- p1 +
    tm_shape(block[block$folds %in% folds_train, ]) +
      tm_borders() + tm_text(text = "folds", palette = "black", size = .7) +
    tm_shape(block[block$folds %in% folds_test, ]) +
      tm_polygons(palette = "grey70", alpha = .7) +
      tm_text(text = "folds", palette = "black", size = .7)

  tmap_arrange(p1, p2, nrow = 1)

}


#' @keywords internal
get_prediction <- function(mod, newdata, type = "response") {
  predict(mod, newdata, type = type, block.size = 1e4, discrete = FALSE) %>%
    as.numeric()
}

#' Spatial Cross-validation of a TGB based PPM and logistic regression
#'
#' @importFrom mgcv bam gam predict.gam
#' @export
cv_sdm_mod <- function(
  fold,
  form,
  add_weights = FALSE,
  select      = TRUE,
  measures    = c("auc", "z_score", "brier"),
  family      = quasi(link = "log", variance = "mu"),
  method      = "fREML",
  discrete    = TRUE,
  mod_fun     = "bam",
  gp_zero     = FALSE) {

  mgcv_args <- list(
    formula  = form,
    family   = family,
    method   = method,
    select   = select,
    discrete = discrete,
    control  = list(trace = TRUE))

  if (gp_zero) {
    test_df <- mutate(test_df, by_gp = 0)
  }

  # calculate model evaluation measures, return data frame with results
  scores <- map(sort(unique(train_df$fold)),
    function(.fold) {

      print(paste0("Iteration ", .fold))
      # newdata for evaluation
      newdata <- filter(train_df, fold == .fold)
      if (gp_zero) {
        newdata <- mutate(newdata, by_gp = 0)
      }
      mgcv_args$data <- filter(train_df, fold != .fold)
      if (add_weights) {
        mgcv_args$weights <- mgcv_args$data$wght
      }
      # fit model
      mod <- do.call(mod_fun, mgcv_args)
      # apply model evaluation measures
      prediction <- get_prediction(mod, newdata)
      truth      <- newdata %>% pull(presence)
      ndf <- data.frame(prediction = prediction, truth = truth)
      res <- map_dbl(measures, ~ do.call(.x, args = list(truth, prediction)))
      names(res) <- measures
      # performance on test data
      prediction_test <- get_prediction(mod, test_df)
      res_test <- map_dbl(measures,
        ~ do.call(.x, args = list(test_df$presence, prediction_test)))
      names(res_test) <- measures
      # return results
      data.frame(
        fold     = .fold,
        measures = rep(measures, 2),
        value    = c(res, res_test),
        type     = rep(c("cv", "test"), each = length(measures)))

    })

  do.call(rbind, scores)

}

remove_duplicates <- function(data, form) {

  keep_vars <- union(all.vars(form),
    c("presence", "merge_year", "fold", "longitude", "latitude"))

  # remove duplidates + missing values
  data %>%
    as.data.frame() %>%
    arrange(merge_year, presence, longitude, latitude) %>%
    distinct(longitude, latitude, merge_year, .keep_all = TRUE) %>%
    select(one_of(keep_vars)) %>%
    na.omit()

}
