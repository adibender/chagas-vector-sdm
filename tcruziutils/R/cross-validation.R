add_blocks <- function(fold, where = "train") {

  polygons <- fold$blocks$blocks
  coords <- sp::coordinates(fold[[where]]) %>% as.data.frame()
  sp::coordinates(coords) <- ~longitude + latitude
  sp::proj4string(coords) <- sp::proj4string(polygons)
  ov <- sp::over(coords, polygons)
  fold[[where]]$layer <- ov$layer

  fold

}

#' Add fold ID to SpatialPointsDataFrame
#'
#' @param An object of class \code{SpatialPointsDataFrame}.
#' @param blocks A SpatialPolygonsDataFrame object. Polygons are spatial
#' blocks created with \code{blockCV}.
#' @param name The name of the column in which the fold ID will be stored.
#'
#' @export
add_folds <- function(
  data,
  blocks,
  name      = "fold",
  spat_form = ~longitude + latitude) {

  coords <- data[all.vars(spat_form)]
  sp::coordinates(coords) <- ~longitude + latitude
  sp::proj4string(coords) <- sp::proj4string(blocks)
  ov <- sp::over(coords, blocks)
  data[[name]] <- as.integer(ov$folds)

  data

}

#' Cross validation for gam models based on folds 1-4
#'
#' @inheritParams get_pmap_ci
#' @inheritParams mgcv::gam
#' @export
cv_gam <- function(
  fold,
  form,
  gamma   = 1L,
  family = binomial(),
  mgcv_args = list(
    select  = TRUE,
    control = list(trace  = FALSE),
    method  = "fREML",
    discrete = TRUE),
  measures = c("auc", "binom_deviance"),
  train_folds = 1:4,
  eval_fold = 5) {

  mgcv_args$form   <- form
  mgcv_args$gamma  <- gamma
  mgcv_args$family <- family

  # remove duplicates (w.r.t. longitude, latitude and year (merge_year))
  train_test <- fold$train %>%
    remove_duplicates(form) %>%
    rbind(fold$test %>% as.data.frame() %>% remove_duplicates(form))
  train_df <- filter(train_test, fold %in% train_folds)
  eval_df <- filter(train_test, fold %in% eval_fold)

  # run CV
  .folds <- unique(train_df$fold) %>% sort()
  message(paste0("gamma: ", mgcv_args$gamma, "\n"))

  res <- map_dfr(
    .folds,
    function(.fold) {

      message(paste0("fold: ", .fold))
      mgcv_args$data <- filter(train_df, fold != .fold) %>%
        mutate(wght = if_else(.data$presence == 1, 1,
          sum(.data$presence) / sum(.data$presence == 0)))
        mgcv_args$weights <- mgcv_args$data$wght
      mod <- try(do.call(mgcv::bam, mgcv_args))
      if (class(mod)[1] == "try-error") {

        candidates <- setdiff(all.vars(mgcv_args$form),
          c("presence", "longitude", "latitude"))
        form <- make_gam_formula(mgcv_args$data[1, ], candidates)
        mgcv_args$formula <- form
        mod <- do.call(mgcv::gam, mgcv_args)

      }
      train_eval <- eval_mod2(mod, newdata = filter(train_df, fold == .fold),
        measures = measures)
      eval <- eval_mod2(mod, newdata = eval_df, measures = measures)
      tibble(
        species = fold$species,
        fold    = .fold,
        train   = list(train_eval),
        test    = list(eval))

    }) %>%
    mutate(gamma = gamma)

}


eval_mod2 <- function(mod, newdata, measures, response_var = "presence") {

  newdata <- add_prediction(mod, newdata, discrete = FALSE)
  res_measures <- map(
    measures,
     ~do.call(.x, args = list(
      newdata[[response_var]],
      newdata[["prediction"]])))
  names(res_measures) <- measures
  res_df <- bind_cols(res_measures) %>%
    mutate(
      eval_df = list(
        tibble(
          response   = newdata[[response_var]],
          prediction = newdata[["prediction"]])))


}


add_formula <- function(settings_df, data) {
  settings_df %>%
  mutate(
    form = map2(type, GP,
      function(.x, .y) {
        form <- make_gam_formula(
          data       = data,
          candidates = candidates,
          f          = "s",
          type       = .x)
        if (.y) {
          form %>% add_gp()
        } else {
          form
        }
      }))
}

get_fold_eval <- function(cv_list, where = "train", what = "auc") {

  cv_list %>%
    bind_rows() %>%
    group_by(gamma) %>%
    mutate(x = map_dbl(.data[[where]], ~.x[[what]])) %>%
    select(gamma, fold, x) %>%
    tidyr::spread(fold, x, sep = "_")

}

#' Cross validation for one setting
#'
#' @export
cv_setting <- function(
  fold,
  type,
  GP,
  candidates,
  env_grids,
  gamma       = 1:4,
  train_folds = 1:4,
  eval_fold   = 5) {

  train_test <- as.data.frame(fold[["train"]]) %>%
    rbind(as.data.frame(fold[["test"]]))
  train <- filter(train_test, .data$fold %in% train_folds)
  test  <- filter(train_test, .data$fold == eval_fold)

  form <- make_gam_formula(train, candidates, type = type)
  if(GP) {
    form <- add_gp(form)
  }
  # if no smooth effects in model, no reason to run model with different penalties
  if(type == "linear" & GP == FALSE) {
    cv_list <- map(gamma[1], ~cv_gam(fold = fold, form = form, gamma = .x,
      train_folds = train_folds, eval_fold = eval_fold))
  } else {
    cv_list <- map(gamma, ~cv_gam(fold = fold, form = form, gamma = .x,
      train_folds = train_folds, eval_fold = eval_fold))
  }

  train_auc_folds <- get_fold_eval(cv_list)
  test_auc_folds  <- get_fold_eval(cv_list, "test")
  opt_gamma <- train_auc_folds %>% tidyr::gather("fold", "auc", -gamma) %>%
    group_by(.data$gamma) %>% summarize(auc = mean(.data$auc)) %>%
    filter(auc == max(.data$auc)) %>% pull(gamma)
  opt_gamma <- opt_gamma[1]

  # refit on all of training data
  train <- train %>%
    mutate(wght = if_else(.data$presence == 1, 1,
      sum(.data$presence) / sum(.data$presence == 0)))
  gam_train <- bam(
    form,
    data     = train,
    family   = binomial(),
    method   = "fREML",
    select   = TRUE,
    discrete = TRUE,
    weights = wght,
    gamma    = opt_gamma)
  smry_gam_train <- summary(gam_train)

  # evaluate on test data
  test_df <- test %>%
    as.data.frame() %>%
    add_prediction(
      mod      = gam_train,
      newdata  = .,
      discrete = FALSE)

  test_auc <- auc(test_df$presence, test_df$prediction)
  # return model and other information
  gam_train$model <- NULL

  tibble(
    type            = type,
    GP              = GP,
    gamma           = opt_gamma,
    test_auc        = test_auc,
    form            = list(form),
    train_auc_folds = list(train_auc_folds),
    test_auc_folds  = list(test_auc_folds),
    mod_smry        = list(smry_gam_train),
    mod             = list(gam_train))

}

#' Cross validation wrapper for multiple settings
#'
#'
#' @export
cv_settings <- function(
  fold,
  settings_df,
  candidates,
  env_grids,
  train_folds = 1:4,
  eval_fold   = 5,
  gamma       = 1:4) {

  res_df <- map_dfr(
    seq_len(nrow(settings_df)),
    function(setting) {
      message(paste0("Setting:\n", settings_df[setting, ]))
      setting <- settings_df[setting, ]
      cv_setting(
        fold        = fold,
        type        = setting$type,
        GP          = setting$GP,
        candidates  = candidates,
        train_folds = train_folds,
        eval_fold   = eval_fold,
        gamma       = gamma)
    })

  # make final prediction on 5x5 raster
  best_ind <- which.max(res_df$test_auc)
  best_mod <- res_df[best_ind,] %>% pull(mod)
  best_mod <- best_mod[[1]]

  ndf  <- env_grids %>% raster::crop(fold$hull) %>%
    raster::mask(fold$hull) %>% grid_to_df()
  ndf  <- ndf %>% add_prediction(mod = best_mod, ., discrete = FALSE)
  grid <- df_to_grid(ndf, env_grids[[1]], "prediction")

  list(species = fold$species, cv_results = res_df, grid = grid)

}

#' smry of cross validation errors
#'
#' @export
eval_smry <- function(cvl, what = "train_auc_folds") {
  cvl[[what]] %>%
    imap_dfr(function(eval_df, index) {
      eval_df %>%
        gather(fold, auc, -gamma) %>%
        group_by(gamma) %>%
        summarize(auc = mean(auc)) %>%
        filter(auc == max(auc)) %>%
        mutate(setting = index)
    })

}

#' Helper functions for cross-validation evaluation
#'
#' @rdname cv_helpers
#' @param cv_list A list of \code{cv_settings} outputs.
#' @param species The species of interest (as character)
#' @export
get_best_mod_smry <- function(cv_list, species)  {

  ind_best <- which.max(cv_list[[species]]$cv_results$test_auc)
  cv_list[[species]]$cv_results$mod_smry[[ind_best]]

}

#' @rdname cv_helpers
#' @inheritParams get_pmap_ci
#' @export
get_best_mod <- function(cv)  {

  ind_best <- which.max(cv$cv_results$test_auc)
  cv$cv_results$mod[[ind_best]]

}

#' @rdname cv_helpers
#' @inheritParams get_pmap_ci
#' @export
get_gp_ls <- function(cv) {

  best_ind <- which.max(cv$cv_results$test_auc)
  mod <- cv$cv_results$mod[best_ind][[1]]
  gp_lgl <- map_lgl(mod$smooth, ~ attr(.x, "class")[1] == "gp.smooth" & .x$by == "NA")
  if (any(gp_lgl)) {
    gp_range <- mod$smooth[[which(gp_lgl)]]$gp.defn[2] * 1e3
  } else {
    gp_range <- NA
  }

  gp_range

}

#' @rdname cv_helpers
#' @inheritParams get_best_mod_smry
#' @param folds A list of objects returned by \code{get_sp_folds}.
#' @param settings_df A data frame of modell settings used when
#' calling \code{cv_settings}.
#' @export
get_spec_smry_df <- function(
  cv_list,
  folds,
  settings_df,
  species) {

  cv       <- cv_list[[species]]
  fold     <- folds[[species]]
  best_ind <- which.max(cv$cv_results$test_auc)
  tibble::tibble(
    species     = species,
    gamma       = cv$cv_results$gamma[best_ind],
    test_auc    = cv$cv_results$test_auc[best_ind],
    gp_range    =  get_gp_ls(cv),
    block_width = fold$block_width,
    range       = fold$vario_range) %>%
  cbind(settings_df[best_ind, ])

}

#' @rdname cv_helpers
#' @inheritParams get_best_mod
#' @export
get_test_auc <- function(cv) {
  round(max(cv$cv_results$test_auc), 2)
}


#' @rdname cv_helpers
#' @inheritParams get_pmap_ci
#' @export
get_best_auc_tab <- function(cv, fold) {

  mod <- get_best_mod(cv)
  test <- fold$test
  pred <- predict(mod, newdata = as.data.frame(test), type = "link",
    discrete = FALSE)
  m_roc <- pROC::roc(response = test$presence, pred, ci = T)

  tibble(
    species = cv$species,
    lower = m_roc$ci[1], mean = m_roc$ci[2], upper = m_roc$ci[3],
    AUC = paste0(
      round(mean, 2),
      " (", round(lower, 2), "; ",
      round(upper, 2), ")"))

}
