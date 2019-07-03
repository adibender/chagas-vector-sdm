# general
library(dplyr)
library(purrr)
# library(pbapply)

# library(tmap)
library(ggplot2)
theme_set(theme_bw())
library(tmap)

# spatial
# library(dismo)
# library(raster)
# devtools::load_all("~/Dropbox/GitHub/mastergrids/")
library(mastergrids)
# library(tcruziutils)
devtools::load_all("../../tcruziutils")

# modeling
library(mgcv)

set.seed(20181205)

## data
# tcruzi_grid <- raster::raster("../../preprocessing/tcruzi_grid.tif")
env_grids <- raster::brick("../../preprocessing/env_grids.grd")
# shp_admin_1 <- readRDS("../../preprocessing/shp_admin_1.Rds")
countries <- readRDS("../../preprocessing/countries.Rds")
p_endemic_zone <- tm_shape(countries) + tm_borders(alpha = .5)
# prediction_df <- grid_to_df(env_grids)


## species data
folds <- readRDS("../../preprocessing/folds_list_vector_presence.Rds")

# set up CV
candidates <- names(env_grids) %>%
  setdiff(c("landcover00", "landcover13", "landcover17", "landcover18"))

# create settings that will be compared
settings_df <- cross_df(list(
  type = c("linear", "smooth", "geospatial"),
  # type = c("linear", "smooth"),
  GP   = c(FALSE, TRUE)))

# f5 <- cv_settings(folds[[5]], settings_df[6 ], candidates, env_grids)

# f5$cv_results$train_auc_folds
# eval_smry(f5$cv_results, "train_auc_folds")
# f5$cv_results$test_auc_folds
# eval_smry(f5$cv_results, "test_auc_folds")
# f5$cv_results$test_auc
# tmap_blocks(folds[[5]])

# p_endemic_zone +
#   tm_shape(folds[[5]]$hull, is.master = TRUE) +
#     tm_polygons(alpha = .1) +
#   tm_shape(f5$grid) +
#     tm_raster(alpha = .8, palette = viridis::magma(1e3),
#     breaks = seq(0, 1, by = .2), style = "cont")


library(parallel)
system.time({
cv_res_all <- mclapply(
  seq_along(folds),
  function(i) {
    species <- folds[[i]]$species
    message(paste0("species: ", species))
    try(cv_settings(fold = folds[[i]], settings_df, candidates, env_grids))
  }, mc.cores = length(folds))
})
names(cv_res_all) <- names(folds)

saveRDS(cv_res_all, "cv-results.Rds")
