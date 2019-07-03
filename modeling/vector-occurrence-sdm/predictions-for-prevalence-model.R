library(purrr)
library(mastergrids)
devtools::load_all("../../tcruziutils")
library(mgcv)
library(tmap)
library(parallel)

tcruzi_grid <- raster::raster("../../preprocessing/tcruzi_grid.tif")
folds <- readRDS("../../preprocessing/folds_list_vector_presence.Rds")
cv_res <- readRDS("cv-results.Rds") %>%
  keep(~class(.x) != "try-error")
folds_iprev <- readRDS("../../preprocessing/folds_infection_prevalence.Rds")
spec_iprev <- folds_iprev$train_test_viprev %>% pull(species) %>% unique()
known_species <- map_chr(cv_res, ~.x$species) %>% intersect(spec_iprev)
cv_res   <- cv_res %>% keep(~.x$species %in% known_species)
folds <- folds %>% keep(~.x$species %in% known_species)

years <- folds_iprev$train_test_viprev %>%
  pull(merge_year) %>%
  unique() %>%
  sort()

dir.create("tgb_predictions")

lapply(years, function(year) {

  common_pred_list <- mclapply(names(cv_res), function(spec) {

    mod <- get_best_mod(cv_res[[spec]])
    cat(paste0("Year: ", year), "\n", paste0("Species: ", spec), "\n\n")
    # load environmental variables for specific year (j)
    env_grids <- raster::brick(paste0("../../preprocessing/env_grids_", year, ".grd"))

    # create a data set with all covariate values from year j for species i
    ndf_year <- env_grids %>%
      grid_to_df(folds[[spec]]$hull, folds[[spec]]$hull)
    ndf_year$species <- spec
    # predicted intensity of thinned PPP (species presence + sampling)
    ndf_year$prediction <- predict(mod, newdata = ndf_year, type = "response",
      # discrete = FALSE,
      block.size = 1e5)
    # prediction as grid/raster
    df_to_grid(ndf_year, tcruzi_grid, column = "prediction")

  }, mc.cores = length(cv_res))
  common_pred_grid <- common_pred_list %>% purrr::reduce(raster::addLayer)
  names(common_pred_grid) <- names(cv_res)
  # save brick with all species for year
  raster::brick(
    common_pred_grid,
    filename = paste0("tgb_predictions/pred_grid_tgb_", year, ".grd"),
    bylayer = FALSE, format = "raster", overwrite = TRUE)

  invisible()

})
