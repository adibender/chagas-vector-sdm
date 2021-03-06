% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross-validation.R
\name{add_folds}
\alias{add_folds}
\title{Add fold ID to SpatialPointsDataFrame}
\usage{
add_folds(data, blocks, name = "fold", spat_form = ~longitude +
  latitude)
}
\arguments{
\item{blocks}{A SpatialPolygonsDataFrame object. Polygons are spatial
blocks created with \code{blockCV}.}

\item{name}{The name of the column in which the fold ID will be stored.}

\item{An}{object of class \code{SpatialPointsDataFrame}.}
}
\description{
Add fold ID to SpatialPointsDataFrame
}
