% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SGATII.R
\name{locationImage}
\alias{locationImage}
\title{Location Density Image}
\usage{
locationImage(
  s,
  xlim,
  ylim,
  nx,
  ny,
  weight = rep_len(1, dim(s)[1L]),
  discard = 0L,
  collapse = TRUE,
  chains = NULL
)
}
\arguments{
\item{s}{a single chain or a list of parallel chains generated by
\code{heliosMetropolis}.}

\item{xlim}{range of the first coordinate.}

\item{ylim}{range of the second coordinate.}

\item{nx}{number of cells in the first coordinate.}

\item{ny}{number of cells in the second coordinate.}

\item{weight}{weights for each location.}

\item{discard}{number of initial samples to discard.}

\item{collapse}{whether to collapse parallel chains to a single chain}

\item{chains}{the set of chains to retain, or \code{NULL}.}
}
\value{
A list with elesments
\item{\code{x}}{the x-ordinates that bound the bins}
\item{\code{y}}{the y-ordinates that bound the bins}
\item{\code{W}}{the weighted image.}
}
\description{
Bin locations to form a 2D density image
}
\details{
Bins the samples for a sequence of locations to produce 2D array
suitable for plotting with \code{image}.  Optionally, a vector of
weights can be provided to differentially weight samples by
location.

This function accepts either a sample from a single mcmc run, or a
list of samples from parallel mcmc runs.  If \code{collapse} is true,
multiple samples are merged and single image is returned,
otherwise an image is returned for each sample.
}
