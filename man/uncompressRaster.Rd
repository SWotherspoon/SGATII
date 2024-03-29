% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SGATII.R
\name{uncompressRaster}
\alias{uncompressRaster}
\title{Uncompress raster vector}
\usage{
uncompressRaster(raster, v)
}
\arguments{
\item{raster}{a template raster}

\item{v}{the vector of finite values}
}
\value{
the vector in raster form.
}
\description{
Uncompress a vector to a raster
}
\details{
The \code{\link[=heliosDiscRasterWeights]{heliosDiscRasterWeights()}} function stores rasters in a
compressed form by recording only the values in the raster that
are finite.  Given a template raster that indicates which values
are finite, this function converts a compressed vector back to
raster form.
}
