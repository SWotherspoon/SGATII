% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SGATII.R
\name{zenithSimulate}
\alias{zenithSimulate}
\title{Solar Zenith Simulation}
\usage{
zenithSimulate(tm, lon, lat, tm.out)
}
\arguments{
\item{tm}{vector of times that specify the template track.}

\item{lon}{vector of longitude that specify the template track.}

\item{lat}{vector of latitude that specify the template track.}

\item{tm.out}{vector of times to which the template is resampled.}
}
\value{
\code{zenithSimulate} returns a data frame with
components
\item{\code{Date}}{times along the simulated track}
\item{\code{Lon}}{longitudes along the simulated track}
\item{\code{Lat}}{latitudes along the simulated track}
\item{\code{Zenith}}{zenith angles along the simulated track}
}
\description{
Simulate zenith angles, along a specified track.
}
\details{
Given times, longitudes and latitudes that specify a template
track, \code{zenithSimulate} interpolates the template onto the new
times specified by \code{tm.out} and computes the solar zenith angle at
each point along the new track.
}
