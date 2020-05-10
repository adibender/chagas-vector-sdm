# Species distribution models of chagas vectors

This folder contains code used to fit species ditribution models for chagas vectors.

**Disclaimer**: The analysis is not reproducible as it relies on
a specific location of covariate layers and some confidential data that
could not be made public. Therefore, we provide here code that demonstrates the
workflow and uses freely available data only. Note: the folder structure is
described further below after the exemplary analysis.

If anything is unclear or not working, please don't hesitate to open up an issue.


# Examplary analysis
This section illustrates the workflow presented in [Bender, Python, Lindsay, et al. (2019)](https://www.biorxiv.org/content/10.1101/738310v1). Note that full details are given in the manuscript and
this code repository. However, as the full data and covariate layers could not be shared, this demo provides a reproducible example with freely available data.
Results based on this data should not differ too much from published results.



```r
# for first run install packages from this repository
# devtools::install("mastergrids", dependencies = TRUE)
# devtools::install("tcruziutils", dependencies = TRUE)
# libraries
devtools::load_all("../mastergrids")
devtools::load_all("../tcruziutils")
library(dplyr)
library(purrr)
# viz
library(ggplot2)
theme_set(theme_bw())
# modeling
library(scam)

# defaults
# this is just the maximum extent of the endemic zone
# used to crop environmental variables, etc.
extent_tcruzi   <- tcruziutils::extent_tcruzi
# country polygons
data(wrld_simpl, package = "maptools")
countries <- wrld_simpl %>% raster::crop(extent_tcruzi)
# colors
Set1   <- RColorBrewer::brewer.pal(9, "Set1")
## set seed as generation of folds is random
set.seed(101850)
```


## Data Import
- In the publication we used the data from
[Ceccarelli, Balsalobre, Medone, et al. (2018)](https://www.nature.com/articles/sdata201871) and some additional data (but not much). The former is openly available ([figshare download link](https://doi.org/10.6084/m9.figshare.c.3946936) (EXCEL file))
- The code below downloads, imports the file and performs some preprocessing



```r
library(httr)
GET("https://ndownloader.figshare.com/files/10302303",
  write_disk(tf <- tempfile(fileext = ".xls")))
```

```
## Response [https://s3-eu-west-1.amazonaws.com/pstorage-npg-968563215/10302303/SciDataData_CitationCeccarellietal.2018.xls]
##   Date: 2020-05-10 16:54
##   Status: 200
##   Content-Type: binary/octet-stream
##   Size: 5.74 MB
## <ON DISK>  /tmp/Rtmpbi1Llj/file32cc611e8598.xls
```

```r
df <- readxl::read_excel(tf, 1L, na = c("", " ", "NR", "NA"))
presence_vector <- df %>%
  mutate(
    area = 25L, # analysis will be performed at 5x5 square km resolution
    Start_year = as.integer(substr(year, 1L, 4L)),
    End_year   = as.integer(substr(year, 6L, 9L)))
presence_vector$reference <- df[[17]]
presence_vector[, 17] <- NULL
presence_vector <- presence_vector %>%
  mutate(
    Public_year = as.integer(stringr::str_extract(.data$reference, "[0-9]{4}")),
    Public_year = ifelse(Public_year > 2018L, NA, Public_year),
    End_year    = ifelse(is.na(End_year), Start_year, End_year)) %>%
  select(scientificName, Start_year, End_year, Public_year, area,
    individualCount, starts_with("decimal"), habitat, reference) %>%
  rename(
    species    = scientificName,
    Latitude   = decimalLatitude,
    Longitude  = decimalLongitude,
    n_observed = individualCount) %>%
  rename_all(~tolower(sub(" ", "_", .))) %>%
  select(-reference)
```

- imputation of the "year" variable if missing



```r
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
```
- remove observations before 2000



```r
presence_vector <- presence_vector %>%
  filter(start_year >= 2000)
```


## Adding spatial/environmental covariates to data set
- in the publication we used raster files available from Servers of
the [Malaria Atlas Project ](https://malariaatlas.org/)
- the environmental variables were extracted for each observation based
on (imputed) year of observation
- here, we use the `raster::getData()` function to obtain the
environmental covariates (we use the same covariate layers for all years)




```r
covs <- raster::getData("worldclim", var = "bio", res = 5)[[1:19]]
env_grids <- raster::crop(covs, extent_tcruzi)
# combine presence + covariate data
covs_ex <- raster::extract(env_grids, as_spatial(presence_vector))
presence_vector <- cbind(presence_vector, covs_ex)
```

## Workflow for one species
The procedure below was iterated over all species (with enough observations) and
consists of the following steps:

1. Create presence/background column for the species of interest
(here `Panstrongylus megistus`)
2. Split data into
  - "training/test" (using spatial blocks/spatial CV), used later for model selection (folds 1-4) and to asses the models "extrapolation performance"
  - "evaluation" (random subsample), used later to asses the models "interpolation performance"
3. The spatial blocks/spatial CV is set up using function
`get_sp_folds` which is a wrapper around the package [**`blockCV`**](https://github.com/rvalavi/blockCV) ([Valavi, Elith, Lahoz-Monfort, et al., 2018](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13107))

4. Fit the model on "training data" (spatial folds 1-4). Here, for illustration, there is no model selection of different models/settings

5. Asses the models performance on data not used during model selection/fit

### 1. Create presence/background dummy


```r
presence_vector <- presence_vector %>%
  mutate(presence = 1L * (species == "Panstrongylus megistus"))
table(presence_vector$presence)
```

```
##
##     0     1
## 13412  1184
```


### 2. Split data


```r
it_pv <- presence_vector %>%
  rsample::initial_split(strata = "species", prop = 4 / 5)
# evaluation data (only used at the very end for "interpolation error" estimate
evaluation_df <- rsample::testing(it_pv)
# train test is split in 5 spatial folds
train_test_df <- rsample::training(it_pv)
```

### 3. Create spatial folds


```r
sp_tt <- as_spatial(train_test_df)
sp_tt <- sp::spTransform(sp_tt, raster::crs(env_grids))
countries <- sp::spTransform(countries, raster::crs(env_grids))
sp_fold_pans_meg <- get_sp_folds(
  data        = sp_tt,
  species     = "Panstrongylus megistus",
  mask        = countries,
  species_var = "species",
  n_blocks    = 50,# number of blocks
  k           = 5, # number of folds
  width       = 5, # width of extended hull (outside of observed )
  calc_range  = FALSE)
```

```
## The best folds was in iteration 84:
##   train_0 train_1 test_0 test_1
## 1    2800     795    485    144
## 2    2574     763    711    176
## 3    2505     706    780    233
## 4    2555     718    730    221
## 5    2706     774    579    165
```


```r
# graphical depiction of spatial folds + presence/background
# within extended hull
tmap_cv(
  sp_fold_pans_meg,
  countries = countries)
```

<p align="center">
<img src="figure/unnamed-chunk-10-1.png" title="plot of chunk unnamed-chunk-10" alt="plot of chunk unnamed-chunk-10" width="400px" style="display: block; margin: auto;" />
</p>

### 4. Fit the model
- In the paper we run a CV for different model specifications on folds 1 - 4
(fit on 3 folds, evaluation on 4th fold)
- best setting/modell w.r.t. average performance for the 4 CV runs is then
refit on folds 1-4 and evaluated using data from fold 5 (this is AUC value
reported in the table)
- this model was also used for evaluation on the random hold-out data (`evaluation_df`) above
- After evaluation the final model was obtained by refitting this model on all
available data (folds 1 -5 and random evaluation data)
- Here we just run one model on folds 1-4, evaluate on fold 5 and hold-out data
and plot the fitted model

#### Fit the GAM

- create `formula` that specifies the linear/additive predictor using
`make_gam_formula`



```r
# create formula
# make_gam_formula: only vars with unique(var) > 20 used
mod_formula <-
  make_gam_formula(
    as.data.frame(sp_fold_pans_meg$train),
    candidates=names(env_grids),
    type = "smooth") %>%
  add_gp() # add gaussian process smooth (see ?mgcv::gp.smooth)
mod_formula
```

```
## presence ~ s(bio1, by = NA) + s(bio2, by = NA) + s(bio3, by = NA) +
##     s(bio4, by = NA) + s(bio5, by = NA) + s(bio6, by = NA) +
##     s(bio7, by = NA) + s(bio8, by = NA) + s(bio9, by = NA) +
##     s(bio10, by = NA) + s(bio11, by = NA) + s(bio12, by = NA) +
##     s(bio13, by = NA) + s(bio14, by = NA) + s(bio15, by = NA) +
##     s(bio16, by = NA) + s(bio17, by = NA) + s(bio18, by = NA) +
##     s(bio19, by = NA) + s(longitude, latitude, bs = "gp")
## <environment: 0x564cd9d23f78>
```


- Fit the GAM using **`mgcv`**
- `mgcv::bam` could be replaced by `mgcv::gam`, but has less demands w.r.t. to
memory reqiurements and offers significant speed-up, especially with `discrete = TRUE` option ([Wood, Li, Shaddick, et al., 2017](https://doi.org/10.1080/01621459.2016.1195744))
- set `discrete = FALSE` when calling the `predict` function to obtain smoother
predictions



```r
# see mgcv::bam for references
mod <- mgcv::bam(
  formula  = mod_formula,
  data     = as.data.frame(sp_fold_pans_meg$train),
  family   = binomial(),
  method   = "fREML", # fast REML
  discrete = TRUE, # speeds up computation
  gamma    = 2L)
```

### 5. Evaluate the model



```r
prediction_test <- predict(
  mod,
  newdata  = as.data.frame(sp_fold_pans_meg$test),
  type     = "response",
  discrete = FALSE)

prediction_eval_df <- predict(
  mod,
  newdata = evaluation_df,
  type = "response",
  discrete = FALSE
)

# evaluation w.r.t. to extrapolation (i.e. fold 5)
MLmetrics::AUC(prediction_test, sp_fold_pans_meg$test$presence)
```

```
## [1] 0.8935416
```


```r
# evaluation w.r.t. interpolation (i.e. random hold-out data)
MLmetrics::AUC(prediction_eval_df, evaluation_df$presence)
```

```
## [1] 0.9632264
```



### 6. Refit model on all data


```r
# extract data points within extended hull of Panstrongylus megistus
df_all <- presence_vector %>%
  as_spatial() %>%
  raster::crop(sp_fold_pans_meg$hull)
# refit model with all data
mod_all <- update(mod, data = as.data.frame(df_all))
```


### 7. Visualize results

- Final prediction:



```r
# newdata with covariate values for each 5x5 pixel within hull of Panstrongylus megistus
env_grids <- env_grids %>%
  raster::crop(sp_fold_pans_meg$hull) %>%
  raster::mask(sp_fold_pans_meg$hull)
ndf <- grid_to_df(env_grids)
# calculate predictions, set discrete = FALSE to obtain smoother predictions
prediction <- predict(mod_all, newdata = ndf, type = "link", discrete = FALSE,
  se = TRUE)
ndf$prediction <- exp(prediction$fit)/(1 + exp(prediction$fit))
# calculate CI
ndf$se <- prediction$se
ci_lower <- prediction$fit - 2*prediction$se
ci_upper <- prediction$fit + 2*prediction$se
ndf$ci_lower <- exp(ci_lower)/(1 + exp(ci_lower))
ndf$ci_upper <- exp(ci_upper)/(1 + exp(ci_upper))
ndf$ci <- ndf$ci_upper - ndf$ci_lower
# retransform df to raster for plotting
pred_raster <- df_to_grid(ndf, env_grids[[1]], "prediction")
tm_shape(raster::crop(countries, raster::extent(sp_fold_pans_meg$hull))) +
  tm_borders() +
  tm_shape(pred_raster) +
  tm_raster(style = "cont", palette = viridis::magma(1e3),
    breaks = seq(0, 1, by = .2), alpha = .8)
```

<p align="center">
<img src="figure/unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" width="400px" style="display: block; margin: auto;" />
</p>

 - Bivariate map (this is not well implemented in the moment in R), manual hacks
 required (alternatively, could predict upper/lower CI and plot CI alongside prediction)




```r
## Note, this is just for illustration. Specific cut-offs and color palette
# for bivariate maps were used in the publication
# create map + legend
# cut points could be specified
bivar_map <- tm_bivariate(ndf, env_grids[[1]], sp_fold_pans_meg)
# draw figure, x and y control position of legend
tm_bivar_draw(bivar_map, x = .55, y = .05)
```

<p align="center">
<img src="figure/unnamed-chunk-17-1.png" title="plot of chunk unnamed-chunk-17" alt="plot of chunk unnamed-chunk-17" width="400px" style="display: block; margin: auto;" />
</p>



# References
NULL


# Session Info

```r
sessionInfo()
```

```
## R version 4.0.0 (2020-04-24)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 18.04.4 LTS
##
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
##
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
##
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base
##
## other attached packages:
##  [1] sf_0.9-3          httr_1.4.1        scam_1.2-5        mgcv_1.8-31
##  [5] nlme_3.1-147      ggplot2_3.3.0     purrr_0.3.4       dplyr_0.8.99.9002
##  [9] tcruziutils_0.0.7 mastergrids_0.0.3 RefManageR_1.2.12 knitr_1.28
##
## loaded via a namespace (and not attached):
##   [1] readxl_1.3.1          spam_2.5-1            backports_1.1.6
##   [4] lwgeom_0.2-3          plyr_1.8.6            lazyeval_0.2.2
##   [7] sp_1.4-1              splines_4.0.0         crosstalk_1.1.0.1
##  [10] listenv_0.8.0         leaflet_2.0.3         gstat_2.0-5
##  [13] usethis_1.6.1         digest_0.6.25         foreach_1.5.0
##  [16] htmltools_0.4.0       viridis_0.5.1         fansi_0.4.1
##  [19] magrittr_1.5          checkmate_2.0.0       memoise_1.1.0
##  [22] MLmetrics_1.1.1       tensor_1.5            remotes_2.1.1
##  [25] recipes_0.1.12        globals_0.12.5        gower_0.2.1
##  [28] xts_0.12-0            rsample_0.0.6         prettyunits_1.1.1
##  [31] colorspace_1.4-1      rgdal_1.4-8           xfun_0.13
##  [34] leafem_0.1.1          callr_3.4.3           crayon_1.3.4
##  [37] jsonlite_1.6.1        spatstat_1.63-3       spatstat.data_1.4-3
##  [40] zoo_1.8-8             survival_3.1-12       iterators_1.0.12
##  [43] glue_1.4.0            stars_0.4-1           polyclip_1.10-0
##  [46] pals_1.6              gtable_0.3.0          ipred_0.9-9
##  [49] pkgbuild_1.0.8        maps_3.3.0            abind_1.4-5
##  [52] scales_1.1.0          mvtnorm_1.1-0         DBI_1.1.0
##  [55] bibtex_0.4.2.2        Rcpp_1.0.4.6          viridisLite_0.3.0
##  [58] progress_1.2.2        units_0.6-6           mapproj_1.2.7
##  [61] dotCall64_1.0-0       Formula_1.2-3         intervals_0.15.2
##  [64] stats4_4.0.0          lava_1.6.7            prodlim_2019.11.13
##  [67] dismo_1.1-4           htmlwidgets_1.5.1     FNN_1.1.3
##  [70] RColorBrewer_1.1-2    geosphere_1.5-10      ellipsis_0.3.0
##  [73] farver_2.0.3          reshape_0.8.8         pkgconfig_2.0.3
##  [76] XML_3.99-0.3          nnet_7.3-14           deldir_0.1-25
##  [79] caret_6.0-86          tidyselect_1.0.0      rlang_0.4.6.9000
##  [82] reshape2_1.4.4        tmaptools_3.0         cellranger_1.1.0
##  [85] munsell_0.5.0         tools_4.0.0           xgboost_1.0.0.2
##  [88] cli_2.0.2             generics_0.0.2        devtools_2.3.0
##  [91] evaluate_0.14         stringr_1.4.0         goftest_1.2-2
##  [94] ModelMetrics_1.2.2.2  processx_3.4.2        leafsync_0.1.0
##  [97] fs_1.4.1              timereg_1.9.4         blockCV_2.1.1
## [100] pec_2019.11.03        pbapply_1.4-2         future_1.17.0
## [103] xml2_1.3.2            compiler_4.0.0        rstudioapi_0.11
## [106] curl_4.3              png_0.1-7             e1071_1.7-3
## [109] testthat_2.3.2        spatstat.utils_1.17-0 spacetime_1.2-3
## [112] tibble_3.0.1          stringi_1.4.6         highr_0.8
## [115] ps_1.3.3              desc_1.2.0            fields_10.3
## [118] rgeos_0.5-3           lattice_0.20-41       Matrix_1.2-18
## [121] classInt_0.4-3        vctrs_0.3.0           pillar_1.4.4
## [124] lifecycle_0.2.0       furrr_0.1.0           pammtools_0.2.2
## [127] cowplot_1.0.0         data.table_1.12.8     raster_3.1-5
## [130] R6_2.4.1              gridExtra_2.3         KernSmooth_2.23-17
## [133] sessioninfo_1.1.1     codetools_0.2-16      dichromat_2.0-0
## [136] MASS_7.3-51.6         assertthat_0.2.1      pkgload_1.0.2
## [139] rprojroot_1.3-2       withr_2.2.0           hms_0.5.3
## [142] parallel_4.0.0        grid_4.0.0            rpart_4.1-15
## [145] timeDate_3043.102     tidyr_1.0.3           class_7.3-17
## [148] automap_1.0-14        tmap_3.0              pROC_1.16.2
## [151] numDeriv_2016.8-1.1   lubridate_1.7.8       base64enc_0.1-3
```



# Folder structure
Overview of project files and folders:

- Preliminary notes:
    + The initial steps of this project requiere access to the **`mastergrids`** folder at Maleria Atlas Project and will thus not be fully reproducible

- **`mastergrids`**:
An **`R`** package that facilitates the import of environmental
raster data from the **`mastergrids`** folder at BDI MAP and some
utility functions to transform rasters to data frames and vice versa
(the data import won't work outside of the BDI/without connection to mastergrids). Also contains two functions `grid_to_df` and `df_to_grid` which
convert RasterLayer/RasterBrick object to a data frame and vice versa.

- **`tcruziutils`**: An **`R`** package that facilitates all steps of the
analysis. Most functions are specific to this project and should not be used
for general purpose projects. It can be loaded at the beginning of the script
using `devtools::load_all("<path to project>/tcruziutils")` or installed via
`devtools::install("<path to project>/tcruziutils")` and then loaded as usual
with `library(tcruziutils)`

- **`infection data`** (not included):
Contains data on infection prevalence and presence in vectors and humans
    + **`External vector database`**: Additional data set (compiled by another
    research group that contains *presence only* data on different
    vector species)

- **`endemic zone`** (not included):
Contains a shape-file that defines the endemic zone of the
disease ("mask"). Can be used to crop environmental grid data and other spatial
objects. The spatial extent defined by this mask is stored in the
**`tcruziutils`** package as `raster::extent` and `sp::bbox` objects for
convenience (see `tcruziutils::extent_tcruzi` and `tcruziutils::bbox_tcruzi`).

- **`polygon boundaries`** (not included):
Shapefiles containing polygon boundaries on administrative district levels 1 and 2. Used to extract location/area and polygon information based on the GAUL code. Cropped according to endemic zone and stored as `shp_admin_1.Rds` and `shp_admin_2.Rds` in the **`preprocessing`** folder

- **`preprocessing`**:
This folder contains the main pre-processing scripts and stores the pre-processed data sets that will be used for modeling

    + `import.R`: Initial data import of the presence/prevalence. For
    observations recorded on a polygon level, also adds the centroid
    and area information of the respective polygon to the data (see function
    `add_spat_dat` in **`tcruziutils`** package)

    + `prep-for-modeling.Rmd`: Builds on the `import.R` script and pre-processes the initially imported data for modeling purposes. This includes:
      - application of inclusion/exclusion criteria
      - imputation of (some) missing data
      - addition of covariate layers to the observed presence/absence data
      (based on coordinates; see `raster::extract` and
      `tcruziutils::add_grid_data`)
      - Splits the data set in *train*/*test* data (and possibly *evaluation*
      data). Additionally creates block-wise cross-validation scheme,
      stored as `fold` variable in the original data.
      - also produces some additional visualizations (stored in
      `tcruzi/figures` as `pdf` and `png`)


- **`modeling`**:

     - **`vector_occurrence-sdm`**:
     Folder containing Species Distribution Models (SDM) based on presense only data:
     - See README therein
