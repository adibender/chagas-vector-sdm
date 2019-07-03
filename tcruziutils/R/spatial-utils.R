#' Transform data frame to spatial object
#'
#' Adds coordinates and projection to data set, turning it into an object
#' of class \code{Spatial*}.
#'
#' @param data Data to transform
#' @importFrom sp coordinates
#' @export
as_spatial <- function(
  data,
  coord_form = ~longitude + latitude,
  proj       = tcruziutils::proj_tcruzi) {

  if (!grepl("Spatial", class(data)[1])) {
    coordinates(data) <- coord_form
  }

  add_projection(data, proj = proj)

}

#' @rdname as_spatial
#' @importFrom raster projection "projection<-"
#' @export
add_projection <- function(
  data,
  proj = tcruziutils::proj_tcruzi) {

  projection(data) <- proj
  data
}

#' Merge data using GAUL code
#'
#' Merges two data sets using GAUL code.
#'
#' @inheritParams sp::merge
#' @param keep Columns from \code{y} that will be kept when merging with \code{x}
#'
#' @importFrom sp merge
#' @export
merge_admin <- function(
  y,
  x,
  keep           = NULL,
  by.x           = "GAUL_CODE",
  by.y           = "admin2_1",
  duplicateGeoms = TRUE,
  all.x          = FALSE) {

  m1 <- sp::merge(x, as.data.frame(y), by.x = by.x, by.y = by.y,
    duplicateGeoms = duplicateGeoms, all.x = all.x)

}



#' Add centroid and area polygon information based on GAUL code.
#'
#' For all \code{candidates} columns in \code{data}, looks up the smallest
#' polygons not \code{NA} (column names have to be of format \code{*x_y},
#' if multiple available per row, will pick the one with largest \code{x} and
#' \code{y}), looks up centroid and area information for respective polygon in
#' \code{shp_admin} data and adds respective information back to \code{data}.
#'
#' @param data A data set containing GAUL codes (possibly on different
#' administrative levels; see \code{candidates}).
#' @param shp_admin A polygon data frame with \code{by} column containing
#' GAUL codes
#' @param by The column in \code{shp_admin} that contains GAUL code
#' @param candidates The columns in \code{data} that contain GAUL codes
#' on different administrative levels.
#'
#' @import dplyr sp
#' @importFrom purrr map map_dbl
#' @importFrom rgeos gCentroid
#' @importFrom raster area
#' @importFrom sp merge split coordinates
#'
#' @export
add_spat_dat <- function(
  data,
  shp_admin,
  keep       = "ID",
  by         = "GAUL_CODE",
  candidates = c("Admin1_1", "Admin2_1"),
  id.var     = "ID") {

  adm <- get_adm_data(data, candidates = candidates)
  adm_shp <- merge(shp_admin, adm, by = "GAUL_CODE", duplicateGeoms = TRUE,
    all.x = FALSE)
  shp_list <- split(adm_shp, f = adm_shp@data[[id.var]])
  centroids_adm <- map(shp_list, ~ gCentroid(.x))
  area_adm <- map_dbl(shp_list, ~ area(.x) / 100000)
  centroids_coords <- do.call(rbind, centroids_adm) %>% coordinates()
  adm <- adm %>%
    select(.data[[id.var]]) %>%
    mutate(
      centroid_x = centroids_coords[, "x"],
      centroid_y = centroids_coords[, "y"],
      area       = area_adm)

  left_join(data, adm)

}


#' @rdname add_spat_dat
#' @inheritParams add_spat_dat
#' @import dplyr
#' @importFrom tidyr gather
#' @keywords internal
get_adm_data <- function(
  data,
  candidates = c("Admin1_1", "Admin2_1"),
  spat_form = ~Longitude + Latitude,
  id_var = "ID") {

  long_var <- all.vars(spat_form)[1]

  admin_cols <- candidates[candidates %in% colnames(data)]
  adm <- data %>%
    filter(is.na(.data[[long_var]])) %>%
    gather(key = "Admin", value = "GAUL_CODE", !!!quo(admin_cols)) %>%
    filter(!is.na(GAUL_CODE))
  adm %>%
    mutate(level1 = substr(Admin, 6, 6)) %>%
    mutate(level2 = substr(Admin, 8, 8)) %>%
    arrange(.data[[id_var]], desc(.data$level1), desc(.data$level2)) %>%
    select(one_of(c(id_var, "Admin", "GAUL_CODE"))) %>%
    group_by(.data[[id_var]]) %>%
    slice(1) %>%
    ungroup()

}



#' Extract GAUL CODE from tcruzi data sets
#'
#'
#' @importFrom tidyr gather
#' @importFrom rlang quo '!!!'
#' @param id_var Name of the variable that can be used as unique identifier.
#' @export
get_gaul_data <- function(data, id_var = "id") {


  admin_cols <- grep("admin._.", colnames(data), value = TRUE)
  adm <- data %>%
    select(one_of(c(id_var, admin_cols))) %>%
    tidyr::gather(key = "admin", value = "GAUL_CODE", !!!(quo(admin_cols))) %>%
    filter(!is.na(GAUL_CODE))

  adm %>%
    mutate(level1 = substr(admin, 6, 6)) %>%
    mutate(level2 = substr(admin, 8, 8)) %>%
    arrange(id, desc(.data$level1), desc(.data$level2)) %>%
    select(one_of(id_var), .data$admin, .data$GAUL_CODE) %>%
    group_by(.data[[id_var]]) %>%
    slice(1) %>%
    ungroup() %>%
    select(one_of(id_var), admin, GAUL_CODE)

}

#' Add column with GAUL CODE to data
#'
#' @importFrom dplyr select left_join
#' @importFrom rlang '!!' enquo
#' @inheritParams as_spatial
#' @inheritParams get_gaul_data
#' @rdname get_gaul_data
#' @export
add_gaul_code <- function(data, id_var = "id") {

  gaul_data <- get_gaul_data(data, id_var = quo_name(id_var)) %>%
    select(-.data$admin)
  data %>% left_join(gaul_data, by = id_var)

}

#' Transform rows with polygon location to multi-rows with grid centroids
#'
#' @importFrom dplyr filter select one_of rename
#' @importFrom purrr map_dfr
#' @importFrom furrr future_map
#' @importFrom raster xyFromCell cellFromPolygon area
#' @importFrom rgeos gCentroid
#' @export
grid_points_in_district <- function(
  data,
  grid,
  admin_polygons,
  gaul_var  = "GAUL_CODE",
  id_var    = "id",
  spat_form = ~longitude + latitude,
  .progress = interactive()) {

  gaul_data <- filter(data, !is.na(.data[[gaul_var]])) %>%
    select(one_of(c(id_var, "GAUL_CODE")))
  gaul_data %>%
    split(f = gaul_data$id) %>%
    furrr::future_map_dfr(function(district) {
      try({
        gaul_code    <- district$GAUL_CODE
        shp_district <- admin_polygons[admin_polygons$GAUL_CODE == gaul_code, ]
        cells <- cellFromPolygon(grid, shp_district)[[1]]
        if (is.null(cells)) {
          if (area(shp_district) / 1e6 < 25) {
            centroids <- gCentroid(shp_district) %>% as.data.frame()
            out <- data.frame(id = district$id, GAUL_CODE = gaul_code) %>%
              cbind(centroids)
          } else {
            out <- data.frame(id = district$id, GAUL_CODE = gaul_code, x = NA,
              y = NA)
          }
        } else {
          out <- cbind.data.frame(
            id        = district$id,
            GAUL_CODE = gaul_code,
            xyFromCell(grid, cells))
        }
        out %>%
          rename(
            !!all.vars(spat_form)[1] := "x",
            !!all.vars(spat_form)[2] := "y")
      })
    },
    .progress = .progress)

}


#' @rdname grid_points_in_district
#' @importFrom rlang quo_name
#' @export
admin_to_multipoint <- function(
  data,
  grid,
  admin_polygons,
  gaul_var = "GAUL_CODE",
  id_var = "id",
  spat_form = ~longitude + latitude,
  .progress = interactive()) {

  multipoint_df <- grid_points_in_district(data, grid, admin_polygons,
    quo_name(gaul_var), quo_name(id_var), spat_form, .progress)
  point_data <- data %>% filter(is.na(.data[[gaul_var]]))
  admin_data <- data %>% filter(!is.na(.data[[gaul_var]])) %>%
    select(-one_of(setdiff(colnames(multipoint_df), id_var))) %>%
    left_join(multipoint_df, by = id_var)
  rbind(point_data, admin_data)

}


#' calculate extent ratio
#'
#'
#' @export
#' @keywords internal
extent_ratio <- function(x, ...) {

  extent_x <- extent(x)
  width    <- extent_x@xmax - extent_x@xmin
  height   <- extent_x@ymax - extent_x@ymin

  height / width

}
