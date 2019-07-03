#' Add quadrature points to presence only data
#'
#' Adds quadrature points (a.k.a. presence background points) to presence only
#' data, possibly stratified by a factor variable(s).
#' In order to do so, covariate information has to be available for all
#' coordinates at which quadrature points are generated (random or
#' structured/gridded). This is usually the case for environmeltal variables
#' available as \code{raster*} objects.
#'
#' @param data The data set that will be augmented with \code{n}
#' quadrature  points.
#' @param queries List of search queries.
#' See \code{\link[mastergrids]{search_files}}.
#' @param hull A polygon object representing the convex hull of the study area.
#' @param year_var The name of the column in which information on year of
#' observation is stored.
#' @param coord_vars The names of the columns in which the longitude and
#' latitude information are stored respectively (in that order).
#' @param strata The name of the column in \code{data} that will be used for
#' stratification. Must be a \code{factor} or \code{character}. After \code{n}
#' quadrature points are sampled, the \code{strata} column will be added and
#' filled by sampling values from the empirical distribution of the
#' \code{strata} variable.
#' @param add_raster_layers Logical. Indicates whether covariate inforamtion should be added in this step.
#' @importFrom raster crop raster
#' @importFrom mastergrids get_grid_data_by_year
#'
#' @export
get_quads <- function(
  data,
  hull,
  queries,
  strata       = NULL,
  year_var     = "start_year",
  coord_vars   = c("longitude", "latitude"),
  add_raster_layers = TRUE,
  ...) {

  quads  <- quad_coords(data, hull, ..., spat_form = as_coord_form(coord_vars))
  quads <- quads %>%
    as.data.frame() %>%
    mutate(id = seq_len(n()) + max(data$id))

  quads <- add_strata(quads, data, strata = c(strata, year_var))
  if (add_raster_layers)  {
    quads <- quads %>%
      get_grid_data_by_year(
      queries,
      year_var  = year_var,
      spat_form = as_coord_form(coord_vars),
      names     = names(queries))
  }

  quads <- quads %>% mutate(presence = 0)
  quads[base::intersect(colnames(data), colnames(quads))]

}

#' @rdname get_quads
#' @param ... Parameters passed to \code{get_quads}
#' @export
add_quads <- function(data, ...) {

  quads <- get_quads(data = data, ...)

  rbind(data[colnames(quads)], quads)

}

#' @keywords internal
#' @importFrom raster xyFromCell values projection "projection<-"
coordinates_from_grid <- function(grid, hull, coord_vars) {

  grid <- crop_grid(grid, hull)
  coords <- raster::xyFromCell(grid, which(!is.na(raster::values(grid)))) %>%
    as.data.frame()
  coordinates(coords) <- ~ x + y
  projection(coords) <- projection(grid)

  coords

}


#' Generate coordinates of background points within a spatial extent
#'
#' @keywords internal
#' @param ... Further arguments passed to \code{\link[spatstat]{quadscheme}}
#' @importFrom raster projection
#' @importFrom spatstat ppp quadscheme default.dummy
#' @inheritParams spatstat::ppp
quad_coords <- function(
  data,
  hull,
  ...,
  spat_form = ~longitude + latitude) {

  lon_var <- all.vars(spat_form)[1]
  lat_var <- all.vars(spat_form)[2]
  win     <- as_owin(hull)
  pp      <- ppp(data[[lon_var]], data[[lat_var]], window = win)
  quad    <- quadscheme(data = pp, ...)
  coords  <- cbind.data.frame(x = quad$dummy$x, y = quad$dummy$y) %>%
    rename(!!lon_var := "x", !!lat_var := "y")
  coords  <- as_spatial(coords, spat_form, proj = projection(hull))
  pol_id <- map_chr(hull@polygons, ~.x@ID)
  hull <- hull[!duplicated(pol_id), ]
  hull_pol <- SpatialPolygons(hull@polygons)
  projection(hull_pol) <- projection(hull)
  ind <- over(coords, hull_pol)

  coords[!is.na(ind), ]

}

subset_ind <- function(ind, method = c("regular", "random"), n = 1e5) {

  method <- match.arg(method)
  n_orig <- length(ind)
  if (n_orig < n) {
    stop("Number of specified sample size larger than number of available
     coordinate sampling points")
  }
  if (method == "regular") {
    sample_ind <- ind[seq(1, n_orig, by = ceiling(max(1, n_orig / n)))]
  } else {
    sample_ind <- sample(ind, n, replace = FALSE)
  }

  sample_ind

}

#' @keywords internal
crop_grid <- function(grid, hull) {
  if (class(grid) != "RasterLayer") {
    return(crop(raster(grid), extent(hull)))
  } else {
    return(crop(grid, extent(hull)))
  }
}

#' @keywords internal
as_coord_form <- function(coord_vars) {
  as.formula(paste0("~", paste0(coord_vars, collapse = "+")))
}

#' @keywords internal
#' @importFrom purrr map map2_dfc
add_strata <- function(quads, data, strata) {
  if (!is.null(strata)) {
    u_strata    <- map(strata, ~ unique(data[[.x]]))
    prob_strata <- map(strata, ~ prop.table(table(data[[.x]])))
    strata_df   <- map2_dfc(u_strata, prob_strata,
      ~ sample(.x, size = nrow(quads), prob = .y, replace = TRUE))
    colnames(strata_df) <- strata
  }
  quads <- cbind(strata_df, quads)
}
