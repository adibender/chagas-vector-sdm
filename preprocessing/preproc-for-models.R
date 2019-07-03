library(dplyr)
library(purrr)
library(furrr)
library(forcats)
# modeling
library(scam)
library(rsample)
# viz
library(ggplot2)
theme_set(theme_bw())
# devtools::install_github("thomasp85/patchwork")
library(patchwork)
library(tmap)
# utils
library(mastergrids)
library(tcruziutils)
data("tcruzi_queries", package = "tcruziutils")

## data sets
# spatial
shp_admin_1 <- readRDS("../preprocessing/shp_admin_1.Rds")
shp_admin_2 <- readRDS("../preprocessing/shp_admin_2.Rds")
countries   <- readRDS("../preprocessing/countries.Rds")

# vector
presence_vector  <- readRDS("../preprocessing/presence_vector.Rds")

## impute missing year dates
presence_vector <- presence_vector %>%
  rename_all(tolower) %>%
  filter(!(is.na(start_year) & is.na(public_year))) %>%
  filter(start_year <= end_year | is.na(start_year) | is.na(end_year)) %>%
  filter(start_year <= public_year | is.na(start_year) | is.na(public_year)) %>%
  filter(start_year >= 2000 | is.na(start_year)) %>%
  filter(public_year >= 2000 | is.na(public_year)) %>%
  filter(end_year <= public_year | is.na(end_year) | is.na(public_year)) %>%
  mutate(
    imputed = is.na(start_year),
    id      = row_number())
# impute start end year
lm_start <- scam(start_year ~  s(public_year, bs = "mpi"),
  data = presence_vector)
lm_end <- scam(end_year ~ s(start_year, bs = "mpi") +
  s(I(public_year - start_year), bs = "mpi"), data = presence_vector)
presence_vector <- presence_vector %>%
  mutate(start_year = ifelse(is.na(start_year),
    as.integer(floor(predict (lm_start, .))), start_year)) %>%
  mutate(end_year = ifelse(is.na(end_year),
    as.integer(ceiling(predict(lm_end, .))), end_year)) %>%
  # few imputed end year larger than public year
  mutate(end_year = pmin(end_year, public_year))  %>%
  # set end_year to start year if still na (not imputed b/c public year = NA)
  mutate(end_year = if_else(is.na(end_year), start_year, end_year))
stopifnot(all(na.omit(presence_vector$start_year <= presence_vector$end_year)))
stopifnot(all(na.omit(presence_vector$end_year <= presence_vector$public_year)))

## add spatial covariate layers (by year if available)
# this won't work outside of BDI/linux/mastergrids mounted at "/media/z/mastergrids"
future::plan("multiprocess", workers = 10)
presence_vector <- presence_vector %>%
  mutate(merge_year = floor(start_year + (end_year - start_year) / 2)) %>%
  as.data.frame() %>%
  get_grid_data_by_year(tcruzi_queries, year_var = "merge_year",
    names = names(tcruzi_queries), .progress = FALSE)
stopifnot(all(presence_vector$start_year <= presence_vector$merge_year))

# create some covariates for later convenience
# - species
# - species2
# - number of species/genus observations
presence_vector <- presence_vector %>%
  mutate(genus = stringr::str_extract(species, "[:alpha:]+")) %>%
  group_by(species) %>%
  mutate(n_species = n()) %>%
  ungroup() %>%
  group_by(genus) %>%
  mutate(n_genus = n()) %>%
  ungroup()
presence_vector <- presence_vector %>%
  mutate(species2 = case_when(
    n_species < 50 ~ paste0(genus, " (other)"),
    TRUE ~ species ))
presence_vector_sp <- presence_vector %>%
  filter(!is.na(longitude)) %>%
  mutate(log_n = log(n_observed))

## training, test and evaluation data sets
set.seed(27112018)
it_pv <- presence_vector %>%
  mutate(presence = 1) %>%
  initial_split(strata = "species2", prop = 4 / 5)
# evaluation data (only used at the very end for the final model)
evaluation_df <- testing(it_pv)
train_test_df <- training(it_pv)
saveRDS(train_test_df, "../preprocessing/presence_vector_train_test.Rds")
saveRDS(evaluation_df, "../preprocessing/presence_vector_evaluation.Rds")


## create spatial CV folds
set.seed(20181203)
future::plan("multiprocess", workers = 10)

tab_species <- table(train_test_df$species2)
# block-wise spatial cross validation for train/test data
included_species <- names(tab_species)[tab_species >= 50]
folds_list <- future_map(included_species,
    function(.x) {
      set.seed(sum(order(stringr::str_split(.x, ""))))
      try(get_sp_folds(
        data        = as_spatial(train_test_df),
        species     = .x,
        mask        = countries,
        rasterLayer = raster::raster("../preprocessing/tcruzi_grid.tif"),
        n_blocks    = 50,
        k           = 5,
        width       = 5))
    })
names(folds_list) <- included_species

folds_list <- folds_list[map_chr(folds_list, ~ class(.x)) != "try-error"]
saveRDS(folds_list, "../preprocessing/folds_list_vector_presence.Rds")
