library(purrr)
library(tmap)
# library(stars)
library(ggplot2)
theme_set(theme_bw())
# library(tcruziutils)
devtools::load_all("../../tcruziutils")
library(dplyr)
library(parallel)

# spatial
# country polygons
countries <- readRDS("../../preprocessing/countries.Rds")
# raster brick with covariate values within endimic zone
env_grids <- raster::brick("../../preprocessing/env_grids.grd")
# island polygons to be excluded from analysis
exclude <- rgdal::readOGR("../../polygon boundaries/islands_americas.shp")

# data held out for evaluation (interpolation performance)
# used here for final model fit (but does not affect evaluation that was based
# on models fit without this data)
ev_df  <- readRDS("../../preprocessing/presence_vector_evaluation.Rds")

# preprocessed data with presence/background data and spatial folds
# each element containt data for one species
folds <- readRDS("../../preprocessing/folds_list_vector_presence.Rds")

# model results ( each element contains results for one species )
# this was used for evaluation
# see cross-val-gam-tgb.Rds
cv_res <- readRDS("cv-results.Rds")

# model specifications used for analysis
settings_df <- purrr::cross_df(list(
  type = c("linear", "smooth", "geospatial"),
  GP   = c(FALSE, TRUE)))

# remove species that errored (and category "other")
ind_keep <- map_lgl(cv_res, ~class(.x)[1] != "try-error") &
  # map_lgl(folds, ~ (sum(.x$train$presence) + sum(.x$test$presence)) >= 100) &
  map_lgl(folds, ~ !grepl("other", .x$species))
cv_res <- cv_res[ind_keep]
folds  <- folds[ind_keep]


## smry table of best settings
settings_best <- map_dfr(
  names(cv_res),
  ~{
    print(.x)
    get_spec_smry_df(cv_res, folds, settings_df, .x)
  })
settings_best <- settings_best %>%
      dplyr::mutate(
        block_gp_ratio    = block_width / gp_range,
        range_gp_ratio    = range / gp_range,
        block_range_ratio = block_width / range)
  settings_best

species_keep <- settings_best %>% pull(species)

cv_res <- cv_res[species_keep]
folds  <- folds[species_keep]


## the best model is refit using all data, including data used previously for
# model validation (fold with blocks numbered 5 and random evaluation data)
pred_list <- mclapply(
  seq_along(cv_res),
  function(ind) {

    .x <- cv_res[[ind]]
    .y <- folds[[ind]]

    print(.y$species)
    hull <- .y$hull
    # mask <- .y$hull %>% raster::crop(exclude)

    ndf <- env_grids %>% grid_to_df(.y$hull, .y$hull)
    ndf$species <- .y$species
    bmod <- get_best_mod(.x)
    form <- formula(bmod)
    ## refit with data from all blocks (1-5)  + random hold-out for final prediction
    # add block 5 to data
    ndata <- rbind(.y$train, .y$test)
    # add hold-out data
    ndata_ev <- as_spatial(ev_df)[.y$hull, ]
    ndata_ev$fold <- NA

    ndata <- rbind(ndata_ev, ndata)
    ndata$presence <- 0
    ndata$presence <- 1L * (ndata$species == .y$species)

    ## readjust weights
    ndata$wght <- 1
    ndata$wght[ndata$presence == 0] <- sum(ndata$presence)/sum(!ndata$presence)
    opt_gamma <- settings_best %>%
      filter(species == .y$species) %>% pull(gamma)
    bmod <- update(bmod, formula = form, data = as.data.frame(ndata),
      gamma = opt_gamma)
    pred <- predict(
      bmod,
      ndf,
      type     = "link",
      discrete = FALSE,
      se       = TRUE)
    ndf$prediction  <- as.numeric(exp(pred$fit) / (1  + exp(pred$fit)))
    ndf$ci_lower    <- as.numeric(exp(pred$fit - 2 * pred$se) /
      (1 + exp(pred$fit - 2 * pred$se)))
    ndf$ci_upper    <- as.numeric(exp(pred$fit + 2 * pred$se) /
      (1 + exp(pred$fit + 2 * pred$se)))
    nas <- rowSums(is.na(ndf[, c("prediction", "ci_lower", "ci_upper")])) > 0
    if( sum(nas) > 0) {
      ndf[nas,] <- NA
    }
    ndf <- ndf %>% filter(!is.na(prediction))
    grid_prediction <- df_to_grid(ndf, env_grids[[1]], column = "prediction") %>%
      raster::mask(exclude, inverse = TRUE)
    grid_cilower    <- df_to_grid(ndf, env_grids[[1]], column = "ci_lower") %>%
      raster::mask(exclude, inverse = TRUE)
    grid_ciupper    <- df_to_grid(ndf, env_grids[[1]], column = "ci_upper") %>%
      raster::mask(exclude, inverse = TRUE)
    names(grid_prediction) <- names(grid_cilower) <- names(grid_ciupper) <-
      .y$species

    gc()

    list(
      species = .y$species,
      final_mod = bmod,
      prediction = grid_prediction,
      ci_lower = grid_cilower,
      ci_upper = grid_ciupper)

  }, mc.cores = 5)

# store the final model within cv_res object
for(i in seq_along(pred_list)) {
  cv_res[[pred_list[[i]]$species]]$final_mod <- pred_list[[i]]$final_mod
}
saveRDS(cv_res, "cv-results-final.Rds")

# combine raster for different species into raster brick
pred_brick <- pred_list %>% map(~.x$prediction) %>% reduce(raster::addLayer)
cil_brick  <- pred_list %>% map(~.x$ci_lower) %>% reduce(raster::addLayer)
ciu_brick  <- pred_list %>% map(~.x$ci_upper) %>% reduce(raster::addLayer)

# these are the final maps available online
raster::brick(
  pred_brick,
  filename  = "predicted-maps-tgb-revision.grd",
  bylayer   = FALSE,
  format    = "raster",
  overwrite = TRUE)
raster::brick(
  cil_brick,
  filename  = "predicted-cil-tgb-revision.grd",
  bylayer   = FALSE,
  format    = "raster",
  overwrite = TRUE)
raster::brick(
  ciu_brick,
  filename  = "predicted-ciu-tgb-revision.grd",
  bylayer   = FALSE,
  format    = "raster",
  overwrite = TRUE)


# save predicted maps (as pdf and png)
for (ind in seq_along(cv_res)) {

    .x <- cv_res[[ind]]
    .y <- folds[[ind]]

    print(.y$species)
    hull <- .y$hull
    spec <- sub(" ", ".", .y$species)
    # mask <- .y$hull %>% raster::crop(exclude)
    shp <- countries %>% raster::crop(raster::extent(hull))
    pred_grid <- raster::brick("predicted-maps-tgb-revision.grd")[[spec]] %>%
      raster::crop(raster::extent(hull))
    pred_map <- tm_shape(shp)  +
      tm_borders(alpha = .8) +
    tm_shape(pred_grid) +
      tm_raster(
        alpha = .7, style = "cont", palette = viridis::magma(1e3),
        breaks = seq(0, 1, by = .2), title = "prediction") +
    tm_layout(
      title             = .y$species,
      title.size        = 2.5,
      title.fontface    = "italic",
      legend.position   = c("left", "bottom"),
      legend.text.size  = 1.2,
      legend.hist.size  = 1.2,
      legend.title.size = 1.5)

    path <- paste0("figures/final-predictions/pred-map-final-",
        sub(" ", "-", .y$species), ".pdf")
    tmap_save(pred_map, filename = path, width = 7, height = 7)

    tmap_save(pred_map,
      filename = paste0("figures/final-predictions/pred-map-final-",
        sub(" ", "-", .y$species), ".png"), width = 7, height = 7)
  }


################################################################################
############################# Bivariate Maps ###################################
################################################################################
# color palette
pal <- matrix(pals::arc.bluepink(), nrow = 4, ncol = 4)
# use some sensible cut offs for prediction
# .5 is baseline because presence/background are reweighted to have equal weights
# thus, probs above indicate presence liklier than absence.
cut_fit = c(0, .25, .5, .75, 1)
# SE of >=.3 implies that lower and upper CI are in different categories of the
# prediction as defined above, then go down in steps of 1/2
cut_se = c(0, .075, .15, .3, 1)

# loop over species and create bivariate maps
for(i in seq_along(cv_res)) {
    cv <- cv_res[[i]]
    species_i <- cv$species
    print(species_i)
    fold <- folds[[species_i]]
    tm_bv <- tm_bivar_raster(
      path_pred = "predicted-maps-tgb-revision.grd",
      path_ciu = "predicted-ciu-tgb-revision.grd",
      path_cil = "predicted-cil-tgb-revision.grd",
      fold,
      palette          = as.vector(t(pal)),
      add_points_after = FALSE,
      cut_fit          = cut_fit,
      cut_se           = cut_se
    )
    tm_bv_drawn <- tm_bivar_draw(tm_bv, x = .75, y = .1, width = .2, height = .2
    )

    path <- paste0("figures/bivariate_maps/map_bivar_", sub(" ", "_", species_i))
    png(paste0(path, ".png"), width = 600, height = 600)
    print(tm_bv_drawn)
    dev.off()

    pdf(paste0(path, ".pdf"), width = 7, height = 7)
    print(tm_bv_drawn)
    dev.off()

  }

### bivar maps for figure 2 with adjusted legend positions
# color palette
pal <- matrix(pals::arc.bluepink(), nrow = 4, ncol = 4)
# use some sensible cut offs for prediction
# .5 is baseline because presence/TGB absence are reweighted to have equal weights
# in sum. probs above indicate presence liklier than absence.
cut_fit = c(0, .25, .5, .75, 1)
# ci width of .3 indicates that upper and lower CI will not be in same category
# of predicted probabilities as defined above, move down by factor 1/2
cut_se = c(0, .075, .15, .3, 1)
spec_figure_2 <- c("Triatoma infestans", "Triatoma dimidiata", "Triatoma gerstaeckeri",
  "Rhodnius pictipes", "Panstrongylus geniculatus")
leg_pos <- list(
  "Triatoma infestans" = c(x = .7, y = .1, width = .3, height = .3),
  "Triatoma dimidiata" = c(x = .2, y = .1, width = .3, height = .3),
  "Triatoma gerstaeckeri" = c(x = .59, y = .2, width = .25, height = .25),
  "Rhodnius pictipes" = c(x = .75, y = .65, width = .2, height = .2),
  "Panstrongylus geniculatus" = c(x = .75, y = .65, width = .25, height = .25)
)

folds <- readRDS("../../preprocessing/folds_list_vector_presence.Rds")[spec_figure_2]

# loop over species for figure 2
for(i in spec_figure_2) {

  fold <- folds[[i]]
  species_i <- fold$species
  print(species_i)

  tm_bv <- tm_bivar_raster(
    path_pred = "predicted-maps-tgb-revision.grd",
    path_ciu = "predicted-ciu-tgb-revision.grd",
    path_cil = "predicted-cil-tgb-revision.grd",
    fold,
    palette          = as.vector(t(pal)),
    add_points_after = FALSE,
    cut_fit          = cut_fit,
    cut_se           = cut_se
  )
  pars_i <- leg_pos[[species_i]]
  tm_bv_drawn <- tm_bivar_draw(tm_bv, x = pars_i["x"], y = pars_i["y"],
    width = pars_i["width"], height = pars_i["height"] )

  path <- paste0("figures/bivariate_maps/map_bivar_", sub(" ", "_", species_i))
  png(paste0(path, ".png"), width = 600, height = 600)
  print(tm_bv_drawn)
  dev.off()

  pdf(paste0(path, ".pdf"), width = 7, height = 7)
  print(tm_bv_drawn)
  dev.off()

}


################################################################################
############################### feature importance #############################
################################################################################
# Supplement figure on feature importance
cv_res <- readRDS("cv-results-final.Rds")
folds <- folds[names(cv_res)]
mean_term_contrib <- purrr::map2(
  cv_res,
  folds,
  ~ {
    print(.x$species)
    ndf <- env_grids %>% grid_to_df(.y$hull, .y$hull)
    terms <- predict(.x$final_mod, newdata = ndf, type = "terms")
    terms <- abs(terms)
    terms <- terms/rowSums(terms)
    # mean_term_contrib <-
    colMeans(terms) * 100
  }
)

nenv <- names(env_grids)

term_contrib <- purrr::map_dfr(
  mean_term_contrib,
  ~{
    out_df <- purrr:::map_dfr(
      nenv,
      function(var) {
        ind <- grepl(var, names(.x))
        if(all(!ind)) {
          val <- NA
        } else {
          val <- .x[ind]
        }
        data.frame(term = var, value = val)
      }
    )
  }, .id = "species")
contrib_gp <- term_contrib %>%
  group_by(species) %>%
  summarize(value = 100 - sum(value,na.rm=T)) %>%
  mutate(term = "GP")
term_contrib <- rbind(term_contrib, contrib_gp) %>%
  mutate(term = factor(
    term,
    levels = c(nenv, "GP"),
    labels = c("accessibility", "elevation", "slope", "water",
      "evergreen needleleaf forest", "evergreen broadleaf forest",
      "deciduous needleleaf forest", "diciduous broadleaf forest", "mixed forest cover",
      "closed shrubland cover", "open shrubland cover", "woody savannah cover",
      "savannah cover", "grassland cover", "wetland cover", "cropland cover",
      "urban and built up area cover", "cropland / vegetation mosaic",
      "snow and ice cover", "barren area cover", "unclassified land cover",
      "no data on land cover", "daytime temperature", "nighttime temperature",
      "diurnal temperature difference", "nighttime lights", "human population count",
      "vegetation index", "wetness on bare soil", "surface wetness",
      "rainfall", "urbanicity", "Gaussian process")))

saveRDS(term_contrib, "mean_term_contrib.Rds")

top3_per_species <- term_contrib %>%
  group_by(species) %>%
  arrange(desc(value), .by_group = TRUE) %>%
  slice(1:3)

readr::write_csv(top3_per_species, "term_contrib_top3_per_species.csv")

heat <- ggplot(term_contrib, aes(x = term, y = species)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colors = viridis::plasma(1e3), limits = c(0, 100)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("figures/heat-map-term-contribution.pdf", width = 10, height = 9)
ggsave("figures/heat-map-term-contribution.png", width = 10, height = 9)

term_contrib_av <- term_contrib %>% group_by(term) %>%
  summarize(mean_contrib = mean(value, na.rm = TRUE)) %>%
  arrange(desc(mean_contrib))

readr::write_csv(term_contrib_av, "term_contrib_av.csv")
