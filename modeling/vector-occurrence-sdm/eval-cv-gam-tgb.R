library(purrr)
devtools::load_all("../../tcruziutils")
library(mgcv)
library(tmap)

folds <- readRDS("../../preprocessing/folds_list_vector_presence.Rds")
cv_res <- readRDS("cv-results.Rds")

# endemic zone country polygons
countries <- readRDS("../../preprocessing/countries.Rds")
p_endemic_zone <- tm_shape(countries) + tm_borders(alpha = .5)

env_grids <- raster::brick("../../preprocessing/env_grids.grd")

p_endemic_zone <- function(countries, hull) {
  countries <- countries %>%
    raster::crop(hull)
}

get_test_auc <- function(cv) {
  round(max(cv$cv_results$test_auc), 2)
}

map_blocks <- function(cv, fold, countries) {

  block       <- fold$blocks$blocks
  folds_train <- fold$train$fold
  folds_test  <- fold$test$fold
  df <- rbind(fold$train, fold$test)
  n_presence  <- sum(df$presence)

  tm_shape(countries) + tm_borders(alpha = .5) +
  tm_shape(fold$hull) + tm_borders(alpha = .5) +
  tm_shape(df[sample(1:nrow(df), nrow(df)), ]) +
    tm_dots(col = "presence", style = "cat", palette = c("steelblue", "black"),
      alpha = .7, size = .15, title = paste0("presence (n = ", n_presence, ")"))  +
  tm_shape(block[block$folds %in% folds_train, ]) +
    tm_borders() + tm_text(text = "folds", palette = "black", size = .7) +
  tm_shape(block[block$folds %in% folds_test, ]) +
    tm_polygons(palette = "grey70", alpha = .7) +
    tm_text(text = "folds", palette = "black", size = .7) +
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

map_sdm <- function(cv, fold, countries) {

  print(cv$species)
  tm_shape(countries) + tm_borders(alpha = .5) +
  tm_shape(fold$hull) + tm_borders(alpha = .5) +
  tm_shape(cv$grid) +
    tm_raster(alpha = .7, breaks = seq(0, 1, by = .2),
      palette = viridis::magma(1e3), style = "cont",
      title = paste0("prediction (AUC: ", get_test_auc(cv), ")")) +
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

tmap_arrange(
  map_blocks(cv_res[[5]], folds[[5]], countries),
  map_sdm(cv_res[[5]], folds[[5]], countries),
  nrow = 1)
ind_keep <- map_lgl(cv_res, ~class(.x)[1] != "try-error") &
  map_lgl(folds, ~ (sum(.x$train$presence) + sum(.x$test$presence)) >= 100) &
  map_lgl(folds, ~ !grepl("other", .x$species))
sdm_map_list <- map2(cv_res[ind_keep], folds[ind_keep],
    ~tmap_arrange(
      map_blocks(.x, .y, countries),
      map_sdm(.x, .y, countries),
      nrow = 2))

pdf("species_maps5.pdf", width = 6, height = 12)
for (map in sdm_map_list) {
  print(map)
}
dev.off()
# individual maps
for (i in seq_along(sdm_map_list)) {
  species <- folds[ind_keep][[i]]$species
  pdf(paste0("prediction-maps/map5x5_", sub(" ", "-", species), ".pdf"), width = 6, height = 12)
  print(map)
  dev.off()
}
for (i in seq_along(sdm_map_list)) {
  print(species)
  species <- folds[ind_keep][[i]]$species
  png(paste0("prediction-maps/map5x5_", sub(" ", "-", species), ".png"), width = 300, height = 600)
  print(sdm_map_list[[i]])
  dev.off()
}

png("species_maps.png", width = 300, height = 600)
for (map in sdm_map_list) {
  print(map)
}
dev.off()

# save predictions as GeoTIFF

dir.create("prediction-maps")
dir.create("presence-data")

for (i in seq_along(cv_res[ind_keep])) {

  print(names(cv_res)[ind_keep][i])

  grid_i <- cv_res[ind_keep][[i]]$grid
  df_i <- subset(
      rbind(
        folds[ind_keep][[i]]$train,
        folds[ind_keep][[i]]$test),
      subset = presence == 1,
      select = c("presence"))
  names(grid_i) <- "prediction"
  raster::writeRaster(
    grid_i,
    filename = paste0("prediction-maps/",
      sub(" ", "-", cv_res[ind_keep][[i]]$species),
      "-prediction.tiff"),
    overwrite = TRUE)
  readr::write_csv(as.data.frame(df_i),
    path = paste0("presence-data/",
      sub(" ", "-", cv_res[ind_keep][[i]]$species), "-presence.csv"))

}

### check ranges and block size
block_widths <- map_dbl(folds, ~.x[["block_width"]])
ranges <- map_dbl(folds, ~.x[["range"]])
sum(block_widths > ranges)
names(folds)[block_widths > ranges]
block_widths / ranges

data.frame(block = block_widths[ind_keep], range = ranges[ind_keep],
  auc = map_dbl(cv_res[ind_keep]))
