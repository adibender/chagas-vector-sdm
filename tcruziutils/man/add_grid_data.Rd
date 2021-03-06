% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid-utils.R
\name{add_grid_data}
\alias{add_grid_data}
\title{Adds environmental variables from grid data to point data}
\usage{
add_grid_data(data, grid, col_name = NULL, year_var = "start_year",
  spat_form = ~Longitude + Latitude)
}
\arguments{
\item{data}{A data frame with spatial coordinates}

\item{grid}{A list of gridded data sets. Length 1 if synoptic data, otherwise
one element for each year for which grid data is available.
Must be a named list. A column with the name of the list will be
added to \code{data}, containing respective grid data information for each
year available in \code{grid}. If \code{data} contains years that are not
available in \code{grid}, the first element of \code{grid} will be used.}

\item{col_name}{The name of the newly created variable can be specified here,
otherwise will be extracted from \code{names(grid)}.}

\item{year_var}{The name of the variable in \code{data} from which the
years will be extracted for which matches will be looked for in \code{grid}.}

\item{spat_form}{Specifies longitude and latitude variables in
\code{data}. See \code{\link[sp]{coordinates}}.}
}
\description{
Adds environmental variables from grid data to point data
}
