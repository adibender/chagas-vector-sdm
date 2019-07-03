library(purrr)
library(tmap)
# library(stars)
library(ggplot2)
theme_set(theme_bw())
# library(tcruziutils)
devtools::load_all("../../tcruziutils")
library(dplyr)

# colors
palette_bin <- c("black", RColorBrewer::brewer.pal(9, "Set1")[2])

# spatial
countries <- readRDS("../../preprocessing/countries.Rds")
env_grids <- raster::brick("../../preprocessing/env_grids.grd")

# data
orig   <- readr::read_csv("../../preprocessing/species-list-unprocessed.csv") %>%
  mutate(genus = stringr::str_extract(species, "[:alpha:]+"))
tt_df  <- readRDS("../../preprocessing/presence_vector_train_test.Rds")
ev_df  <- readRDS("../../preprocessing/presence_vector_evaluation.Rds")
tot_df <- rbind(tt_df, ev_df)
# total number of data points used for analysis
nrow(tot_df)
testthat::expect_identical(nrow(tot_df), 15215L)
# unique species
tab_species <- tot_df %>% group_by(species) %>% summarize(n = n()) %>%
  arrange(desc(n))
# unique genera
tab_genus <- tot_df %>% group_by(genus) %>% summarize(n = n()) %>%
  arrange(desc(n))

folds <- readRDS("../../preprocessing/folds_list_vector_presence.Rds")

# model results
cv_res <- readRDS("cv-results.Rds")

# model specifications used for analysis
settings_df <- cross_df(list(
  type = c("linear", "smooth", "geospatial"),
  # type = c("linear", "smooth"),
  GP   = c(FALSE, TRUE)))

# summary table, number of presence/background points in train/test data
n_df <- map_dfr(folds,
  ~tibble(
    species     = .x$species,
    n_train     = sum(.x$train$presence),
    n_train_tgb = sum(.x$train$presence == 0),
    n_test      = sum(.x$test$presence),
    n_test_tgb  = sum(.x$test$presence == 0),
    n_total     = n_train + n_test,
    n_total_tgb = n_train_tgb + n_test_tgb))

# remove species that errored (and category "other")
ind_keep <- map_lgl(cv_res, ~class(.x)[1] != "try-error") &
  # map_lgl(folds, ~ (sum(.x$train$presence) + sum(.x$test$presence)) >= 100) &
  map_lgl(folds, ~ !grepl("other", .x$species))
cv_res <- cv_res[ind_keep]
folds <- folds[ind_keep]


## smry table of best settings
settings_best <- map_dfr(names(cv_res),
    ~get_spec_smry_df(cv_res, folds, settings_df, .x)) %>%
  dplyr::mutate(
    block_gp_ratio    = block_width / gp_range,
    range_gp_ratio    = range / gp_range,
    block_range_ratio = block_width / range)
settings_best
# species for which block width was smaller than autocorrelation range
filter(settings_best, block_range_ratio < 1)

# species for which block width was smaller than GP range
filter(settings_best, block_gp_ratio < 1)

species_keep <- settings_best %>% pull(species)

cv_res <- cv_res[species_keep]
folds  <- folds[species_keep]


#### save brick with predicted values
pred_list <- mclapply(
  seq_along(cv_res),
  function(ind) {
    .x <- cv_res[[ind]]
    .y <- folds[[ind]]
    ndf <- env_grids %>% grid_to_df(.y$hull, .y$hull)
    ndf$species <- .y$species
    pred <- predict(
      get_best_mod(.x),
      ndf,
      type     = "link",
      discrete = FALSE,
      se       = TRUE)
    ndf$prediction <- exp(pred$fit) / (1  + exp(pred$fit))
    ndf$ci_lower   <- exp(pred$fit - 2 * pred$se) /
      (1 + exp(pred$fit - 2 * pred$se))
    ndf$ci_upper   <- exp(pred$fit + 2 * pred$se) /
      (1 + exp(pred$fit + 2 * pred$se))
    grid_prediction <- df_to_grid(ndf, env_grids[[1]], column = "prediction")
    grid_cilower    <- df_to_grid(ndf, env_grids[[1]], column = "ci_lower")
    grid_ciupper    <- df_to_grid(ndf, env_grids[[1]], column = "ci_upper")
    names(grid_prediction) <- names(grid_cilower) <- names(grid_ciupper) <-
      .y$species

    list(
      prediction = grid_prediction,
      ci_lower = grid_cilower,
      ci_upper = grid_ciupper)
}, mc.cores = lenght(cv_res))

pred_brick <- pred_list %>% map(~.x$prediction) %>% reduce(raster::addLayer)
cil_brick  <- pred_list %>% map(~.x$ci_lower) %>% reduce(raster::addLayer)
ciu_brick  <- pred_list %>% map(~.x$ci_upper) %>% reduce(raster::addLayer)
raster::brick(
  pred_brick,
  filename  = "predicted-maps-tgb.grd",
  bylayer   = FALSE,
  format    = "raster",
  overwrite = TRUE)
raster::brick(
  cil_brick,
  filename  = "predicted-cil-tgb.grd",
  bylayer   = FALSE,
  format    = "raster",
  overwrite = TRUE)
raster::brick(
  ciu_brick,
  filename  = "predicted-ciu-tgb.grd",
  bylayer   = FALSE,
  format    = "raster",
  overwrite = TRUE)


#### Figures
# species of interest
species_main <- c(
  "Panstrongylus geniculatus", "Panstrongylus megistus",
  "Triatoma barberi", "Triatoma brasiliensis", "Triatoma infestans",
  "Triatoma pseudomaculata")
species_br_ratio <- filter(settings_best, block_range_ratio < 1) %>%
  pull(species)
species_rest <- setdiff(species_keep, species_main)

## spatial blocking
p_spatial_blocks <- map(species_keep,
  ~ {
    print(.x)
    tmap_blocks(folds[[.x]], countries, palette = rev(palette_bin))
  })
names(p_spatial_blocks) <- species_keep

pdf("figures/spatial-blocks.pdf", width = 8, height = 4)
for(block in p_spatial_blocks) {
  print(block)
}
dev.off()

## Figure illustrating constuction of background points and CV setup
fold        <- folds[["Panstrongylus megistus"]]
blocks      <- fold$blocks$blocks
folds_train <- unique(fold$train$fold)
folds_test  <- unique(fold$test$fold)

data_presence <- rbind(
  fold$train[fold$train$presence == 1, ],
  fold$test[fold$test$presence == 1, ])
data_all <- rbind(fold$train, fold$test)
sn    <- sample(1:nrow(data_all), nrow(data_all), replace = FALSE)
hull  <- rgeos::gConvexHull(data_presence)
ehull <- rgeos::gBuffer(hull, width = 5)
hull  <- hull %>%
  raster::intersect(countries)
ehull <- ehull %>%
  raster::intersect(countries)
p_initial_hull <- tm_shape(countries) +
  tm_borders(alpha = .1) +
  tm_shape(hull, is.master = TRUE) +
    tm_polygons(alpha = .4, border.alpha = .1) +
  tm_shape(data_presence) + tm_dots() +
  tm_layout(title = "(A)")
p_extended_hull <- tm_shape(countries) +
  tm_borders(alpha = .1) +
  tm_shape(hull) +
    tm_polygons(alpha = .2, border.alpha = .1) +
  tm_shape(ehull, is.master = TRUE) +
    tm_polygons(alpha = .2, border.alpha = .1) +
  tm_shape(data_presence) + tm_dots() +
  tm_layout(title = "(B)")
p_final <- tm_shape(countries) +
  tm_borders(alpha = .1) +
  tm_shape(hull) +
    tm_polygons(alpha = .2, border.alpha = .1) +
  tm_shape(ehull, is.master = TRUE) +
    tm_polygons(alpha = .2, border.alpha = .1) +
  tm_shape(data_all[sn, ]) +
    tm_dots(col = "presence", style = "cat", legend.show = FALSE,
    palette = rev(c("black", "steelblue")), alpha = .7) +
  tm_add_legend("symbol", col = c("steelblue", "black"), size = c(1, 1),
    labels = c("background", "presence")) +
  tm_layout(title = "(C)", legend.text.size = c(1.2))
p_final2 <-
  tm_shape(countries) +
    tm_borders(alpha = .2) +
  tm_shape(hull) +
    tm_polygons(alpha = .2, border.alpha = .1) +
  tm_shape(ehull, is.master = TRUE) +
    tm_polygons(alpha = .2, border.alpha = .1) +
  tm_shape(data_all[sn, ]) +
    tm_dots(col = "presence", style = "cat", legend.show = FALSE,
    palette = rev(c("black", "steelblue")), alpha = .5) +
  tm_shape(blocks[blocks$folds %in% folds_train, ]) +
    tm_borders(alpha = .7) +
  tm_text(text = "folds", palette = "black", size = 0.8, fontface = "bold") +
  tm_shape(blocks[blocks$folds %in% folds_test, ]) +
    tm_polygons(col="grey70", alpha = .7) +
  tm_text(text = "folds", palette = "black", size = .8, fontface="bold") +
  tm_add_legend("symbol", col = c("steelblue", "black"), size = c(1, 1),
    labels = c("background", "presence")) +
  tm_layout(title = "(D)", legend.text.size = c(1.2))

p_block_generation <- tmap_arrange(p_initial_hull, p_extended_hull, p_final,
  p_final2, nrow = 1, ncol = 4)
tmap_save(p_block_generation, "figures/p-block-generation.png",
  width = 9, height = 3, dpi = 900)
tmap_save(p_block_generation, "figures/p-block-generation.pdf",
  width = 9, height = 3)

## bivariate maps -> prediction + uncertainty
bivar_maps_main <- map(species_main,
  ~tm_bivar_cv(cv_res[[.x]], folds[[.x]], env_grids, palette = tolochko.redblue(9)))
legend_pos <- list(
  c(.1, .1, .25, .25),
  c(.5, .1, .25, .25),
  c(.15, .15, .25, .25),
  c(.2, .1, .25, .25),
  c(.15, .15, .25, .25),
  c(.5, .1, .25, .25))
bivar_plots <- map2(bivar_maps_main, legend_pos,
  ~ tm_bivar_draw(.x, .y[1], .y[2], .y[3], .y[4]))
pdf("figures/bivar_maps_main.pdf", width = 9, height = 12)
gridExtra::grid.arrange(grobs = bivar_plots, nrow = 3, ncol = 2)
dev.off()

# predictive maps + CI
# returns tmap objects containing prediction, lower ci and upper ci maps
library(profvis)
pr <- profvis({get_pmap_ci(cv_res[[2]], folds[[2]], env_grids, countries,
  width = 10)})
ci_maps <- map2(cv_res, folds,
  ~{
    print(.x$species)
    get_pmap_ci(.x, .y, env_grids, countries, width = 10)
  })
names(ci_maps) <- names(cv_res)

tmaps_pred <- map(ci_maps, ~.x$pred)
names(tmaps_pred) <- names(ci_maps)

p_pred_1 <- tmap_arrange(tmaps_pred[species_main], ncol = 2, nrow = 3,
  outer.margins = 0)

tmap_save(p_pred_1, "figures/p_pred_main.pdf", width = 9, height = 12)
tmap_save(p_pred_1, "figures/p_pred_main.png", width = 9, height = 12, dpi = 800)
# tmap_save(p_pred_2, "figures/p_pred_2.png", width = 9, height = 12, dpi = 800)
# tmap_save(p_pred_3, "figures/p_pred_3.png", width = 9, height = 12, dpi = 800)
# tmap_save(p_pred_4, "figures/p_pred_4.png", width = 9, height = 12, dpi = 800)
# tmap_save(p_pred_5, "figures/p_pred_5.png", width = 9, height = 12, dpi = 800)

pdf("figures/maps_with_ci.pdf", width = 12, height = 4)
for(ci_map in ci_maps) {
  print(tmap_arrange(ci_map, nrow = 1, outer.margins = 0))
}
dev.off()

#### Bivariate maps for prediction and uncertainty
bivar_maps_main <- map(
  species_main[1],
  ~tm_bivar_cv(cv_res[[.x]], folds[[.x]], env_grids,
    palette = pals::tolochko.redblue(9)))

library(ggpubr)
ggarrange(bivar_maps_main[[1]][[1]], bivar_maps_main[[2]][[1]],
  nrow = 1L, ncol = 2, align = "hv")

#### additional evaluations
# Summary table for publication
## number train/test
n_df <- n_df %>% filter(species %in% species_keep)
n_df %>% print(n = 50)

## AUC table
auc_tab <- map_dfr(
    cv_res,
    ~tibble(species = .x$species, AUC = get_test_auc(.x))) %>%
  left_join(n_df) %>%
  left_join(select(settings_best, "species", "type", "GP", "gamma")) %>%
  mutate(ms =
    case_when(
      type == "linear" ~ 1L,
      type == "smooth" ~ 3L,
      TRUE             ~ 5L) + 1L*(GP)) %>%
  mutate(
    gamma = as.character(gamma),
    gamma =
      case_when(
        ms == 1 ~ "-",
        TRUE ~ gamma)) %>%
  arrange(desc(n_total))

summary(auc_tab$AUC)
sd(auc_tab$AUC)

## auc tab based on random hold out data
auc_oos <- readRDS("auc-oos-eval-df.Rds")
auc_tab <- left_join(auc_tab, auc_oos, by = "species")

auc_tab_print <- auc_tab %>%
  select(species, n_total, n_test, n_total_tgb, n_test_tgb, ms, gamma, AUC.x, AUC.y) %>%
  mutate(
    species = case_when(
      species %in% species_br_ratio ~ paste0(species, "*"),
      TRUE ~ species ),
    species = paste0("\\emph{", species, "}"),
    n_total = paste0(n_total, " (", n_test, ")"),
    n_total_tgb = paste0(n_total_tgb, " (", n_test_tgb, ")"),
    ms = paste0(ms, " (", gamma, ")"),
    AUC = paste0(AUC.x, " (", round(AUC.y, 2), ")")) %>%
  select(-n_test, -n_test_tgb, -gamma, -AUC.x, -AUC.y) %>%
  rename("n (presence)" = n_total, "n (background)" = n_total_tgb,
    "Model specification" = ms)
## latex table for publication
auc_tab_print %>%
  xtable::xtable() %>%
  print(include.rownames = FALSE)

auc_tab_main <- auc_tab %>% filter(species %in% species_main)
auc_tab_main %>%
  select(species, n_total, n_total_tgb, AUC) %>%
  xtable::xtable() %>%
  print(include.rownames = FALSE)
auc_tab_rest <- auc_tab %>% filter(species %in% species_rest)

auc_tab_rest %>%
  select(species, n_total, n_total_tgb, AUC) %>%
  xtable::xtable() %>%
  print(include.rownames = FALSE)


# auc_tab2 <- map2_dfr(cv_res, folds, ~get_best_auc_tab(.x, .y))
