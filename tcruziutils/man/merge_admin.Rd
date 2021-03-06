% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatial-utils.R
\name{merge_admin}
\alias{merge_admin}
\title{Merge data using GAUL code}
\usage{
merge_admin(y, x, keep = NULL, by.x = "GAUL_CODE", by.y = "admin2_1",
  duplicateGeoms = TRUE, all.x = FALSE)
}
\arguments{
\item{y}{object of class \code{data.frame}, or any other class that can be coerced to a data.frame with \code{as.data.frame}}

\item{x}{object deriving from \link{Spatial}}

\item{keep}{Columns from \code{y} that will be kept when merging with \code{x}}

\item{by.x}{specifications of the common columns. See 'Details' in (base) \code{\link[base]{merge}}.}

\item{by.y}{specifications of the common columns. See 'Details' in (base) \code{\link[base]{merge}}.}

\item{duplicateGeoms}{logical; if \code{TRUE} geometries in \code{x} are
  duplicated if there are multiple matches between records in \code{x}
  and \code{y}}

\item{all.x}{logical; if TRUE, then the returned object will have all rows of x, even those that has no matching row in y. These rows will have NAs in those columns that are usually filled with values from y}
}
\description{
Merges two data sets using GAUL code.
}
