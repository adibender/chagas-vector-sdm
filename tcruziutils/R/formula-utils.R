#' Create gam formula depending on number of unique values in covariates
#'
#' @importFrom purrr map_int imap_dfr map_chr
#' @importFrom dplyr tibble
#' @param data The data set that will be modeled. Needs to have a \code{fold}
#' column indicating cross-validation folds.
#' @param candidates Which variables from \code{data} to consider when
#' building the model formula.
#' @param lhs The left hand side of the formula.
#' @param bs The type of basis functions to use when
#' \code{type = "smooth"} or \code{type = "geospatial"}
#' @param f The name of the smooth term constructor.
#' @param Names of longitude and latitude variables in \code{data}.
#' @param type The type of linear predictor to construct. One of \code{linear},
#' \code{smooth} or \code{geospatial}.
#' @param n_min Only variables that have at least \code{n_min} unique values
#' in each cross-validation fold will be included in the formula.
#' @export
#' @keywords internal
make_gam_formula <- function(
  data,
  candidates    = NULL,
  lhs           = "presence ~",
  bs            = "tp",
  f             = "s",
  long_lat_vars = c("longitude", "latitude"),
  type          = c("linear", "smooth", "geospatial"),
  n_min         = 20L,
  by = NA) {

  type <- match.arg(type)
  candidates       <- preproc_candidates(data, candidates)
  candidates       <- setdiff(candidates, long_lat_vars)
  candidate_groups <- classify_groups(data, candidates) %>% filter(n >= n_min)

  terms <- candidate_groups %>% pull(var)
  if (type == "smooth") {
    terms <- terms %>% map_chr(~paste0(f, "(", .x, ", by = ", by, ")"))
  }
  if(type == "geospatial") {
    terms <- terms %>%
      map_chr(~paste0(f, "(", long_lat_vars[1], ", ", long_lat_vars[2],
        ", by = ", .x, ", bs = 'gp')"))
  }
  if (length(terms) > 0) {
    lhs %>% paste(paste0(terms, collapse = " + ")) %>% as.formula()
  } else {
    lhs %>% paste0(1) %>% as.formula()
  }

}


#' Add Gaussian process term to model formula
#'
#' @rdname make_gam_formula
add_gp <- function(formula) {

  update(as.formula(formula), .~. + s(longitude, latitude, bs = "gp"))

}


preproc_candidates <- function(data, candidates = NULL) {

  if (is.null(candidates)) {
    colnames(data)
  } else {
    intersect(candidates, colnames(data))
  }

}

classify_groups <- function(data, candidates) {

  unique_df <- map_dfr(candidates,
    ~ tibble(n = min_unique(data, var = .x), var = .x))
  unique_df %>%
    mutate(type =
      case_when(
        n < 10 ~ "linear",
        TRUE ~ "smooth"))

}


min_unique <- function(data, var) {

  data %>%
    group_by(fold) %>%
    summarize(n = length(unique(.data[[var]]))) %>%
    pull(n) %>%
    min()

}

#https://stevencarlislewalker.wordpress.com/2012/08/06/merging-combining-adding-together-two-formula-objects-in-r/
#' @export
merge.formula <- function(form1, form2, ...) {

  # get character strings of the names for the responses
  # (i.e. left hand sides, lhs)
  lhs1 <- deparse(form1[[2]])
  lhs2 <- deparse(form2[[2]])
  if(lhs1 != lhs2) stop('both formulas must have the same response')

  # get character strings of the right hand sides
  rhs1 <- attr(terms(form1), "term.labels")
  rhs2 <- attr(terms(form2), "term.labels")

  # create the merged rhs and lhs in character string form
  rhs <- c(rhs1, rhs2)
  lhs <- lhs1

  # put the two sides together with the amazing
  # reformulate function
  out <- reformulate(rhs, lhs)

  # set the environment of the formula (i.e. where should
  # R look for variables when data aren't specified?)
  environment(out) <- parent.frame()

  return(out)
}
