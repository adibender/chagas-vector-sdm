#' Predict model on grid points within hull, set values outside to 0
#'
#' @importFrom mastergrids grid_to_df df_to_grid
#' @importFrom mgcv predict.gam predict.bam
#' @inheritParams add_probs
#' @param fold A fold object, containing species information
#' @param mod A model object (model fit to data in \code{fold}.
#' @export
ppm_to_prob_raster <- function(fold, mod, pred_grid) {

  hull <- fold$hull
  pred_grid_sub <- grid_to_df(pred_grid, mask = hull)
  pred_grid_sub$prediction <- predict(mod$model, pred_grid_sub, type = "response")
  # pred_grid_sub$prob <- pred_grid_sub$prediction / sum(pred_grid_sub$prediction)
  df_to_grid(pred_grid_sub, pred_grid[[1]], "prediction", outside = 0)

}

#' Create multinomialy weighted observations from raster cells within polygon
#'
#' @importFrom raster cellFromPolygon xyFromCell extract
#' @importFrom purrr map_dfr
#' @importFrom dplyr mutate
#' @param multipoint_data An augmented data set where polygon-level observations
#' are represented by multiple points (centroids of 5x5 grid cells within
#' the polygon)
#' @param folds A folds object (list) that contains fold information for
#' each species
#' @param mods A list of models fit to respective \code{folds}.
#' @param pred_grid The grid to which \code{models} will be predicted.
#' @export
add_probs <- function(multipoint_data, folds, mods, pred_grid) {

  prob_rasters <- furrr::future_map2(folds_pres, mods,
  ~ppm_to_prob_raster(.x, .y, env_grids), .progress = interactive())


  prob_raster <- ppm_to_prob_raster(fold, mod, pred_grid)
  admin_codes <- data@data[, c("admin1_1", "admin2_1")] %>%
    mutate(GAUL_CODE = ifelse(!is.na(admin2_1), admin2_1, admin1_1))
  prob_df <- map_dfr(admin_codes$GAUL_CODE, function(district) {
    shp <- shp_admin[shp_admin$GAUL_CODE == district, ]
    cells <- cellFromPolygon(prob_raster, shp)[[1]]
    coords_district <- xyFromCell(prob_raster, cells)
    out <- cbind.data.frame(coords_district,
      prob = extract(prob_raster, coords_district)) %>%
    mutate(
      prob_multinom = prob / sum(prob),
      GAUL_CODE     = district)

  })


}
