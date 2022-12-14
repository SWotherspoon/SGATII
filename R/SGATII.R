##' These functions wrap and unwrap a sequence of longitudes around
##' the dateline.
##'
##' The `wrapLon` function wraps the longitudes back into the interval
##' [lmin,lmin+360).  The `unwrapLon` function unwraps a sequence of
##' longitudes so the the initial point lies in [lmin,lmin+360), but
##' the subsequent longitudes in the sequence may wander outside that
##' range.
##'
##' @title Wrap Locations Around the Dateline.
##' @param lon a vector of longitudes
##' @param lmin western boundary for wrapped longitudes
##' @return a vector of longitudes
##' @export
wrapLon <- function(lon,lmin=-180)
  (lon-lmin)%%360+lmin

##' @rdname wrapLon
##' @export
unwrapLon <- function(lon,lmin=-180)
  cumsum(c(wrapLon(lon[1L],lmin),wrapLon(diff(lon))))



## Solar Zenith/Sunrise/Sunset calculations
##
## The functions presented here are based on code and the excel
## spreadsheet from the NOAA site
##
##       http://www.esrl.noaa.gov/gmd/grad/solcalc/
##


##' Calculate solar time, the equation of time and solar declination
##'
##' The solar time, the equation of time and the sine and cosine of
##' the solar declination are calculated for the times specified by
##' `tm` using the same methods as
##' <www.esrl.noaa.gov/gmd/grad/solcalc/>.
##'
##' @title Solar Time and Declination
##' @param tm a vector of POSIXct times.
##' @return A list containing the following vectors.
##' \item{`solarTime`}{the solar time (degrees)}
##' \item{`eqnTime`}{the equation of time (minutes of time)}
##' \item{`sinSolarDec`}{sine of the solar declination}
##' \item{`cosSolarDec`}{cosine of the solar declination}
##' @seealso [zenith()]
##' @examples
##' ## Current solar time
##' solar(Sys.time())
##' @export
solar <- function(tm) {

  rad <- pi/180

  ## Time as Julian day (R form)
  Jd <- as.numeric(tm)/86400.0+2440587.5

  ## Time as Julian century [G]
  Jc <- (Jd-2451545)/36525

  ## The geometric mean sun longitude (degrees) [I]
  L0 <- (280.46646+Jc*(36000.76983+0.0003032*Jc))%%360


  ## Geometric mean anomaly for the sun (degrees) [J]
  M <- 357.52911+Jc*(35999.05029-0.0001537*Jc)

  ## The eccentricity of earth's orbit [K]
  e <- 0.016708634-Jc*(0.000042037+0.0000001267*Jc)

  ## Equation of centre for the sun (degrees) [L]
  eqctr <- sin(rad*M)*(1.914602-Jc*(0.004817+0.000014*Jc))+
    sin(rad*2*M)*(0.019993-0.000101*Jc)+
      sin(rad*3*M)*0.000289

  ## The true longitude of the sun (degrees) [M]
  lambda0 <- L0 + eqctr

  ## The apparent longitude of the sun (degrees) [P]
  omega <- 125.04-1934.136*Jc
  lambda <- lambda0-0.00569-0.00478*sin(rad*omega)


  ## The mean obliquity of the ecliptic (degrees) [Q]
  seconds <- 21.448-Jc*(46.815+Jc*(0.00059-Jc*(0.001813)))
  obliq0 <- 23+(26+(seconds/60))/60

  ## The corrected obliquity of the ecliptic (degrees) [R]
  omega <- 125.04-1934.136*Jc
  obliq <- obliq0 + 0.00256*cos(rad*omega)

  ## The equation of time (minutes of time) [U,V]
  y <- tan(rad*obliq/2)^2
  eqnTime <- 4/rad*(y*sin(rad*2*L0) -
                      2*e*sin(rad*M) +
                      4*e*y*sin(rad*M)*cos(rad*2*L0) -
                      0.5*y^2*sin(rad*4*L0) -
                      1.25*e^2*sin(rad*2*M))

  ## The sun's declination (radians) [T]
  solarDec <- asin(sin(rad*obliq)*sin(rad*lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)

  ## Solar time unadjusted for longitude (degrees) [AB!!]
  ## Am missing a mod 360 here, but is only used within cosine.
  solarTime <- ((Jd-0.5)%%1*1440+eqnTime)/4
  #solarTime <- ((Jd-2440587.5)*1440+eqnTime)/4

  ## Return solar constants
  list(solarTime=solarTime,
       eqnTime=eqnTime,
       sinSolarDec=sinSolarDec,
       cosSolarDec=cosSolarDec)
}



##' Calculate the solar zenith angle for given times and locations
##'
##' `zenith` uses the solar time and declination calculated by `solar`
##' to compute the solar zenith angle for given times and locations,
##' using the same methods as <www.esrl.noaa.gov/gmd/grad/solcalc/>.
##' This function does not adjust for atmospheric refraction see
##' [refracted()].
##'
##' `cosZenith` is an optimized version of `zenith` that calculates
##' the cosine of the zenith angle.
##' @title Solar Zenith Angle
##' @param sun list of solar time and declination computed by `solar`.
##' @param lon vector of longitudes.
##' @param lat vector latitudes.
##' @return A vector of solar zenith angles (degrees) for the given
##'   locations and times.
##' @seealso [solar()]
##' @examples
##' ## Approx location of Sydney Harbour Bridge
##' lon <- 151.211
##' lat <- -33.852
##' ## Solar zenith angle for noon on the first of May 2000
##' ## at the Sydney Harbour Bridge
##' s <- solar(as.POSIXct("2000-05-01 12:00:00","EST"))
##' zenith(s,lon,lat)
##' ## Cosine of the zenith angle
##' cosZenith(s,lon,lat)
##' @export
zenith <- function(sun,lon,lat) {

  rad <- pi/180

  ## Suns hour angle (degrees) [AC!!]
  hourAngle <- sun$solarTime+lon-180
  #hourAngle <- sun$solarTime%%360+lon-180

  ## Cosine of sun's zenith [AD]
  cosZenith <- (sin(rad*lat)*sun$sinSolarDec+
                cos(rad*lat)*sun$cosSolarDec*cos(rad*hourAngle))

  ## Limit to [-1,1] [!!]
  cosZenith[cosZenith > 1] <- 1
  cosZenith[cosZenith < -1] <- -1

  ## Ignore refraction correction
  acos(cosZenith)/rad
}


##' @rdname zenith
##' @export
cosZenith <- function(sun,lon,lat) {

  rad <- pi/180

  ## Suns hour angle (degrees) [AC!!]
  hourAngle <- sun$solarTime+lon-180
  #hourAngle <- sun$solarTime%%360+lon-180

  ## Cosine of sun's zenith [AD]
  sn <- sin(rad*lat)
  (sn*sun$sinSolarDec+sqrt(1-sn^2)*sun$cosSolarDec*cos(rad*hourAngle))
}



##' Adjust the solar zenith angle for atmospheric refraction.
##'
##' Given a vector of solar zeniths computed by [zenith()],
##' `refracted` calculates the solar zeniths adjusted for the effect
##' of atmospheric refraction.
##'
##' `unrefracted` is the inverse of `refracted`. Given a (single)
##' solar zenith adjusted for the effect of atmospheric refraction,
##' `unrefracted` calculates the solar zenith as computed by
##' [zenith()].
##'
##' @title Atmospheric Refraction
##' @param zenith zenith angle (degrees) to adjust.
##' @return vector of zenith angles (degrees) adjusted for atmospheric
##'   refraction.
##' @examples
##' ## Refraction causes the sun to appears higher on the horizon
##' refracted(85:92)
##' ## unrefracted gives unadjusted zenith
##' unrefracted(refracted(90))
##' @export
refracted <- function(zenith) {
  rad <- pi/180
  elev <- 90-zenith
  te <- tan((rad)*elev)
  ## Atmospheric Refraction [AF]
  r <- ifelse(elev>85,0,
              ifelse(elev>5,58.1/te-0.07/te^3+0.000086/te^5,
                     ifelse(elev>-0.575,
                            1735+elev*(-518.2+elev*(103.4+elev*(-12.79+elev*0.711))),
                            -20.772/te)))
  ## Corrected Zenith [90-AG]
  zenith-r/3600
}




##' @rdname refracted
##' @importFrom stats uniroot
##' @export
unrefracted <- function(zenith)
  uniroot(function(x) refracted(x)-zenith,c(zenith,zenith+2))



##' Distances along a path
##'
##' Compute the great circle distances (in km) between successive
##' locations along path.
##'
##' @title Distance along a path
##' @param x a two column matrix of (lon,lat) locations along the path.
##' @return vector of interpoint distances (km)
##' @export
trackDist <- function(x) {
  n <- nrow(x)
  rad <- pi/180
  cosx2 <- cos(rad*x[,2L])
  sinx2 <- sin(rad*x[,2L])

  6378.137*acos(pmin.int(
    cosx2[-n]*cosx2[-1L]*cos(rad*(x[-1L,1L]-x[-n,1L]))+sinx2[-n]*sinx2[-1L],1))
}


##' Bearing changes along a track
##'
##' Compute the change in bearing between successive locations along
##' path.
##'
##' @title Distance along a path
##' @param x a two column matrix of (lon,lat) locations along the path.
##' @return vector of changes in bearing (degrees)
##' @export
trackBearingChange <- function(x) {
  n <- nrow(x)
  rad <- pi/180
  cosx2 <- cos(rad*x[,2L])
  sinx2 <- sin(rad*x[,2L])


  ## Bearing from one x to the next
  bs <- atan2(sin(rad*(x[-1L,1L]-x[-n,1L]))*cosx2[-1L],
              cosx2[-n]*sinx2[-1L]-
                sinx2[-n]*cosx2[-1L]*cos(rad*(x[-1L,1L]-x[-n,1L])))/rad
  ## Difference bs and fold difference into [-180,180)
  wrapLon(bs[1L-n]-bs[-1L])
}



##' Simulate zenith angles, along a specified track.
##'
##' Given times, longitudes and latitudes that specify a template
##' track, `zenithSimulate` interpolates the template onto the new
##' times specified by `tm.out` and computes the solar zenith angle at
##' each point along the new track.
##'
##' @title Solar Zenith Simulation
##' @param tm vector of times that specify the template track.
##' @param lon vector of longitude that specify the template track.
##' @param lat vector of latitude that specify the template track.
##' @param tm.out vector of times to which the template is resampled.
##' @return `zenithSimulate` returns a data frame with
##' components
##' \item{`Date`}{times along the simulated track}
##' \item{`Lon`}{longitudes along the simulated track}
##' \item{`Lat`}{latitudes along the simulated track}
##' \item{`Zenith`}{zenith angles along the simulated track}
##' @importFrom stats approx
##' @export
zenithSimulate <- function(tm,lon,lat,tm.out) {
  ## unwrap longitudes
  lon <- unwrapLon(lon)
  ## Interpolate track
  keep <- !is.na(lon)
  lon.out <- approx(tm[keep],lon[keep],tm.out,rule=2)$y
  keep <- !is.na(lat)
  lat.out <- approx(tm[keep],lat[keep],tm.out,rule=2)$y
  ## Compute zenith angles
  z <- zenith(solar(tm.out),lon.out,lat.out)
  data.frame(Date=tm.out,
             Lon=lon.out,
             Lat=lat.out,
             Zenith=z)
}



##' Helios Model Structures.
##'
##' Helios requires a model structure that describes the model being
##' fitted by providing functions that compute the contributions to
##' the log posterior from each of the segments along the track
##' (`logp.s`) and each of the locations on the track (`logp.p`),
##' together with a vector of `time` times for which locations are
##' estimated, a two column matrix `x0` of initial estimates of those
##' locations and a logical vector `fixedx` that indicates which of
##' those locations are fixed. This function generate basic model
##' structures that should provide a suitable starting point for most
##' analyses.
##'
##' The `heliosModel` function constructs a model structure that
##' combines the `heliosLightModel` with a behavioural model that
##' assume the average speed of travel between successive locations is
##' Gamma distributed.
##'
##' The helios light model assumes the light sensor cannot detect
##' light when it is dark, but may fail to detect light when it is
##' light because the sensor is obscured.  The model assumes the
##' length of time the sensor can be obscured is exponentially
##' distributed, with rate `alpha`.
##'
##' The `forbid` argument represents the likelihood of light being
##' observed when it is dark.  When `forbid=-Inf` the initial
##' locations `x0` must be consistent with the observed light - there
##' can be no light observed when it is dark on the implied track.
##' Setting `forbid` to a large negative value relaxes this
##' constraint, making the observation of light in the dark unlikely
##' but not impossible.
##'
##' The behavioural model assumes the speed of travel is Gamma
##' distributed.  The shape and rate of the Gamma distribution of
##' speeds is specified by the `beta` parameter. If `beta` is a two
##' element vector the same shape and rate are applied across all
##' segments, but if `beta` is a two column matrix the shape and rate
##' can be specified on a segment by segment basis. By default, the
##' speed of travel is calculated based on the time intervals between
##' the estimated locations (in hours), but the intervals of time
##' actually available for travel can be specified directly with the
##' `dt` argument.  The behavioural model makes a contribution to the
##' posterior on each segment of track.
##'
##' Both the light model and the behavioural model and make a
##' contribution to the posterior on each segment of track.  The user
##' is free to provide a function `logp.s0` of a single argument `x`
##' (the locations) that returns an additional contribution to the
##' posterior for each track segment, and a funtion `logp.p0` that
##' returns an additional contribution to the posterior for each
##' location.
##'
##' Some tags record the maximum light level observed in the
##' preceeding sampling interval, so that when the sun is going down,
##' the tag reports the light from the beginning of the sampling
##' interval. The `max.adjust` argument tries to compensate for this
##' behaviour, and for some tags will yield more reliable estimates.
##' However, this optiion does not impact the along track zenith
##' angles calculated by `predict`, and for some applications these
##' may need manual adjustment.
##'
##' @title Helios Model Structures
##' @param date POSIXct vector of times for day/night observations
##' @param light logical vector of day/night observations
##' @param time times at which to fit locations
##' @param x0 initial estimates of fitted locations.
##' @param alpha rate parameter for sensor obscuration.
##' @param beta parameters of the behavioural model.
##' @param logp.p0 function to evaluate any additional contribution to
##'   the log posterior from the estimated locations on the track
##' @param logp.s0 function to evaluate any additional contribution to
##'   the log posterior from the estimated segments of track
##' @param fixedx logical vector indicating which estimated locations
##'   to hold fixed.
##' @param dt (optional) vector of time intervals (hours) for speed
##'   calculation.
##' @param zenith the solar zenith angle that defines day/night.
##' @param forbid the log likelihood for forbidden light observations
##' @param max.adjust adjust for tags that record the maximum light
##'   level observed in the interval 
##' @return The `heliosModel` function returns list with components
##' \item{`time`}{the times at which locations are estimated}
##' \item{`x0`}{a two column array of initial location estimates.}
##' \item{`fixedx`}{a logical vector indicating which locations are
##'   fixed.}
##' \item{`logp.s`}{function to evaluate the contributions to the log
##'   posterior from the track segments}
##' \item{`logp.p`}{function to evaluate the contributions to the log
##'   posterior from the track points}
##' \item{`logp.l`}{function to evaluate the contributions to the log
##'   posterior from the light model on each segment}
##' \item{`logp.b`}{function to evaluate the contributions to the log
##'   posterior from the behavioural model on each segment}
##' \item{`logp.l0`}{an optimized form of `logp.l` that assumes the
##'   track is stationary.}
##' \item{`predict`}{function to predict the zenith angle along a given
##'   track.}
##' \item{`Date`}{POSIXct vector of times for day/night observations.}
##' \item{`Light`}{logical vector of day/night observations.}
##' The `heliosLightModel` function returns just the `predict` and
##' `logp.l` components.
##' @importFrom stats dgamma
##' @export
heliosModel <- function(date,light,time,x0,
                        alpha,beta,
                        logp.p0=function(x) rep.int(0L,nrow(x)),
                        logp.s0=function(x) rep.int(0L,nrow(x)-1L),
                        fixedx=FALSE,dt=NULL,zenith=96,forbid=-Inf,
                        max.adjust=FALSE) {

  ## Times (hours) between locations
  if(is.null(dt))
    dt <- diff(as.numeric(time)/3600)

  ## Ensure beta is always a matrix
  if(!is.matrix(beta)) beta <- t(beta)

  ## Trim data to location estimates
  light <- light[date >= min(time) & date <= max(time)]
  date <- date[date >= min(time) & date <= max(time)]

  ## Check fixed points
  if(any(is.na(x0[fixedx,]))) warning("Fixed points cannot be NA")

  ## Set fixed segments
  fixeds <- if(length(fixedx)==1L) FALSE else fixedx[-1L] & fixedx[-length(fixedx)]
  
  ## Contribution to log posterior for the light data
  model <- heliosLightModel(date,light,time,alpha,zenith,forbid,max.adjust)
  logp.l <- model$logp.l
  logp.l0 <- model$logp.l0
  predict <- model$predict

  ## The behavioural model
  logp.b <- function(x) {
    spd <- pmax.int(trackDist(x), 1e-06)/dt
    dgamma(spd,beta[,1L],beta[,2L],log=TRUE)
  }

  ## Total contribution on the segments
  logp.s <- function(x) {
    r <- logp.s0(x)+logp.l(x)+logp.b(x)
    r[fixeds] <- 0
    r
  }

  ## Contribution from the locations
  logp.p <- function(x) {
    r <- logp.p0(x)
    r[fixedx] <- 0
    r
  }

  ## Create model
  model <- list(
    ## Location estimates
    time=time,
    x0=x0,
    fixedx=fixedx,
    ## Components of the posterior
    logp.s=logp.s,
    logp.p=logp.p,
    logp.l=logp.l,
    logp.b=logp.b,
    logp.l0=logp.l0,
    ## light prediction
    predict=predict,
    ## Trimmed data
    date=date,
    light=light,
    forbid=forbid
  )
  model
}

##' @rdname heliosModel
##' @export
heliosLightModel <- function(date,light,time,alpha,zenith,forbid,max.adjust) {

  ## Convert times to solar time
  s <- solar(date)

  ## Precalculate cos(z)
  cosZ <- cos(pi/180*zenith)

  segments <- as.numeric(cut(date,time,labels=FALSE,include.lowest=TRUE))
  ks <- c(which(diff(segments)==1),length(date))
  ls <- light > 0

  predict <- function(x) {
    ## Interpolate
    lon <- approx(time,x[,1L],date,rule=2)$y
    lat <- approx(time,x[,2L],date,rule=2)$y
    ## Calculate cosz
    cosZ <- cosZenith(s,lon,lat)
    ## Limit to [-1,1] [!!]
    cosZ[cosZ > 1] <- 1
    cosZ[cosZ < -1] <- -1
    
    ## Ignore refraction correction
    180*acos(cosZ)/pi
  }

  if(max.adjust) {
    
    logp.l <- function(x) {
      ## Interpolate
      lon <- approx(time,x[,1L],date,rule=2)$y
      lat <- approx(time,x[,2L],date,rule=2)$y
      ## Calculate where should be day
      day <- cosZenith(s,lon,lat) > cosZ
      day <- c(FALSE,day[-length(day)]) | day 
      ## Calculate costs per segment
      ifelse(diff(c(0,cumsum(!day & light)[ks]))>0,forbid,-alpha*diff(c(0,cumsum(day & !light)[ks])))
    }
    
    logp.l0 <- function(x) {
      ## Calculate where should be day
      day <- cosZenith(s,x[1L],x[2L]) > cosZ
      day <- c(FALSE,day[-length(day)]) | day 
      ## Calculate costs per segment
      ifelse(diff(c(0,cumsum(!day & light)[ks]))>0,forbid,-alpha*diff(c(0,cumsum(day & !light)[ks])))
    }
    
    list(predict=predict,logp.l=logp.l,logp.l0=logp.l0)

  } else {
    
    logp.l <- function(x) {
      ## Interpolate
      lon <- approx(time,x[,1L],date,rule=2)$y
      lat <- approx(time,x[,2L],date,rule=2)$y
      ## Calculate where should be day
      day <- cosZenith(s,lon,lat) > cosZ
      ## Calculate costs per segment
      ifelse(diff(c(0,cumsum(!day & light)[ks]))>0,forbid,-alpha*diff(c(0,cumsum(day & !light)[ks])))
    }
    
    logp.l0 <- function(x) {
      ## Calculate where should be day
      day <- cosZenith(s,x[1L],x[2L]) > cosZ
      ## Calculate costs per segment
      ifelse(diff(c(0,cumsum(!day & light)[ks]))>0,forbid,-alpha*diff(c(0,cumsum(day & !light)[ks])))
    }
    
    list(predict=predict,logp.l=logp.l,logp.l0=logp.l0)
  }

    
}


##' Metropolis sampler for Helios
##'
##' These functions draw samples from the posterior for the Helios
##' model by the Metropolis algorithm.  If `parallel` is `TRUE` the
##' chains are drawn in parallel with [future_lapply()].
##'
##' The starting values `x0` should be specified as a two column
##' matrix, or as a list of two column matrices, one for each chain.
##'
##' @title Metropolis Samplers
##' @param model a model structure as generated by [`heliosModel()`]
##'   or equivalent.
##' @param proposal function for drawing proposals for x.
##' @param x0 starting values for estimated locations x.
##' @param iters number of samples to draw.
##' @param thin rate at which to thin samples.
##' @param chains number of chains to sample.
##' @param parallel enable parallel evaluation through futures.
##' @return If there are r samples drawn for each of q chains of p
##'   parameters at n locations, Helios will return a list containing
##'   \item{`model`}{the model structure}
##'   \item{`x`}{a list of n x p x r arrays of estimated locations from
##'              the q chains}
##' @importFrom stats runif
##' @importFrom future.apply future_lapply
##' @export
heliosMetropolis <- function(model,proposal,x0=NULL,
                              iters=1000L,thin=10L,chains=1L,
                              parallel=TRUE) {

  ## Initialize x0
  if(is.null(x0)) x0 <- model$x0
  ## Expand starting values for multiple chains
  x0 <- rep(if(is.list(x0)) x0 else list(x0),length.out=chains)

  ## Number of locations
  n <- nrow(x0[[1L]])
  ## Number of parameters
  m <- ncol(x0[[1L]])

  ## Extract model components
  logp.s <- model$logp.s
  logp.p <- model$logp.p
  fixedx <- model$fixedx


  ## PARALLEL - parallelise this loop
  sampler <- function(k1) {
    ## Allocate chain
    ch <- array(0,c(n,m,iters))

    ## Initialize
    x1 <- x0[[k1]]
    ## Drop dimnames for speed
    dimnames(x1) <- NULL

    ## Contribution to logp from each track segment, padded with
    ## fictitious segments at each end
    logp.s1 <- c(0,logp.s(x1),0)

    ## Contributions to the log posterior from the locations
    logp.p1 <- logp.p(x1)

    for(k2 in 1L:iters) {
      for(k3 in 1L:thin) {

        ## Propose all x at once, keeping fixed points
        xp <- proposal(x1)
        xp[fixedx,] <- x1[fixedx,]

        for(rb in 1L:2L) {

          ## Modify every second location
          is <- seq.int(rb,n,by=2L)
          x2 <- x1
          x2[is,] <- xp[is,]

          ## Contributions to the log posterior from each track
          ## segment, padded with a fictious segment at each end.
          logp.s2 <- c(0,logp.s(x2),0)

          ## Contributions to the log posterior from the locations
          logp.p2 <- logp.p(x2)

          ## The contribution to the log posterior from each location are
          ## the sums of the contributions from the surrounding segments.
          logp1 <- logp.s1[is]+logp.s1[is+1L]+logp.p1[is]
          logp2 <- logp.s2[is]+logp.s2[is+1L]+logp.p2[is]

          ## MH rule - compute indices of the accepted points.
          accept <- is[logp2-logp1 > log(runif(length(is)))]
          x1[accept,] <- x2[accept,]
          logp.p1[accept] <- logp.p2[accept]
          logp.s1[accept] <- logp.s2[accept]
          logp.s1[accept+1L] <- logp.s2[accept+1L]
        }
      }
      ## Store the current state
      ch[,,k2] <- x1
    }
    ch
  }
  if(parallel)
    chains <- future_lapply(seq_len(chains),sampler,future.seed=TRUE)
  else
    chains <- lapply(seq_len(chains),sampler)
  list(model=model,x=chains)
}


##' Construct the Helios weights on a discrete spatial grid for use
##' with 'heliosDiscMinimalTrack' and `heliosDiscSampler`.
##'
##' The [`heliosDiscMinimalTrack()`] and [`heliosDiscSampler()`] functions
##' generate approximate tracks on a discrete spatial grid that can be
##' used to initialize [`heliosMetropolis()`] sampler.  Both these
##' functions require a weights object that defines the spatial grid
##' and the relative likelihood of the points on this grid based on
##' the observed light data.
##'
##' To generate a weights object with `heliosDiscGridWeights` the user
##' must specify the grid from which locations will be drawn as a two
##' column matrix of (lon,lat) points. The grid can be arbitrary but
##' must be cover the range of the track.
##'
##' The `heliosDiscRasterWeights` function provides an alternative
##' interface in which the grid is specified as a terra raster.  If
##' that raster has values, locations with values that not finite are
##' excluded from the grid.  The `extractHeliosDiscRasterWeights`
##' function extracts the weights for a single track point as a
##' raster, and can be used to visualize the computed weights.
##'
##' The weights are stored as a matrix with a row for each grid point
##' and column for each point along the track, and the entries of each
##' column are the relative posterior probabilities that that track
##' point passes through the corresponding grid point.
##'
##' If `points=FALSE` the weights are calculated from the light data
##' along tracks segments and the contribution to the posterior from
##' the track points.  If `points=FALSE` the weights are calculated
##' from the contribution to the posterior from the light data alone.
##' If the contribution to the posterior from the tracks points is
##' only used to exclude locations, it is more efficient to exclude
##' those points from the grid directly and set `points=FALSE`.
##'
##' If `autoscale` is `NA` the scaling of the weights is determined
##' directly by the posterior. Otherwise the weights are scaled so the
##' maximum weight is 1 and the minimum weight of a feasible grid
##' point is the value of `autoscale`.  If scaling is a vector, its
##' elements define the scaling used at each individual location.
##'
##' @title Helios Discrete Grid Weights
##' @param grid a two column matrix of (lon,lat) grid points.
##' @param model a model structure as generated by [`heliosModel()`]
##'   or equivalent.
##' @param autoscale a lower bound for the (finite) weights
##' @param points should the contribution to the posterior from the
##'   track points be included.
##' @param raster a terra raster used to define the grid.
##' @param k the rack point ot extract
##' @param log show the weights be returned on the log scale.
##' @param weights a weights object generated with
##'   `heliosDiscGridWeights` or `heliosDiscRasterWeights`
##' @return The `heliosDiscGridWeights` and `heliosDiscRasterWeights`
##'   functions return a list with the grid of locations, the model
##'   object and the matrix of weights. The
##'   `extractHeliosDiscRasterWeights` function returns a column of
##'   weights as a raster.
##' @export
heliosDiscGridWeights <- function(grid,model,
                                  autoscale=NA,points=TRUE) {

  m <- nrow(grid)
  n <- length(model$time)
  forbid <- model$forbid
  autoscale <- rep(autoscale,length.out=n)

  ## Matrix of contributions to the posterior for each point of the
  ## grid, assuming animal is stationary.
  W <- matrix(0,m,n)
  if(points) {
    ## Include contribution from points
    for(k in seq_len(m)) {
      lp <- model$logp.l0(grid[k,])
      W[k,] <- c(0,lp) + c(lp,0) + model$logp.p(grid[rep(k,n),])
    }
  } else {
    ## Assume forbidden points are excluded from the grid
    for(k in seq_len(m)) {
      lp <- model$logp.l0(grid[k,])
      W[k,] <- c(0,lp) + c(lp,0)
    }
  }

  ## Standardize at each location
  for(k in seq_len(n)) {
    w <- W[,k]-max(W[,k])
    if(!is.na(autoscale[k])) {
      mn <- min(w[is.finite(w) & w > forbid/2])
      if(mn < 0) w <- (log(autoscale[k])/mn)*w
    }
    W[,k] <- w
  }
  W <- exp(W)

  list(grid=grid,model=model,W=W)
}




##' @rdname heliosDiscGridWeights
##' @importFrom terra hasValues values xyFromCell ncell
##' @export
heliosDiscRasterWeights <- function(raster,model,
                                    autoscale=NA,points=FALSE) {


  if(hasValues(raster))
    grid <- xyFromCell(raster,which(is.finite(values(raster))))
  else
    grid <- xyFromCell(raster,seq_len(ncell(raster)))

  heliosDiscGridWeights(grid,model,autoscale,points)
}


##' @rdname heliosDiscGridWeights
##' @importFrom terra deepcopy hasValues values values<-
##' @export
extractHeliosDiscRasterWeights <- function(raster,weights,k,log=FALSE) {

  r <- deepcopy(raster)
  w <- weights$W[,k]
  if(log) w <- log(w)
  if(hasValues(r))
    values(r)[is.finite(values(r))] <- w
  else
    values(r) <- weights$W[,k] <- w
  r
}



##' Draw samples for Helios from a discrete spatial grid.
##'
##' The `heliosDiscSampler` function samples from a Helios model where
##' the locations are constrained to a discrete spatial grid, and can
##' be used to generate an initial track for the
##' [`heliosMetropolis()`] sampler.  The sampler requires a weights
##' object generated with either [`heliosDiscGridWeights()`] or
##' [`heliosDiscRasterWeights()`] with `log=FALSE`.
##'
##' The weights define the proposal distribution used by the sampler
##' and represent the relative probabilities of a grid point being
##' selected by the proposal distribution at location along the track.
##' The weights are derived from the model likelihood assuming the tag
##' is stationary.
##'
##' This sampler does not require initial location estimates and it
##' only uses the initial locations `x0` from the model to specify
##' fixed points, the remaining elements of `x0` are ignored.
##'
##' This is an independence sampler in which the proposal distribution
##' is based on an approximation the to helios likelihood obtained by
##' assuming the tag is stationary. So for very fast moving
##' individuals, the method can perform poorly.
##'
##' @title Discrete Sampler
##' @param weights a weights object generated with
##'   `heliosDiscGridWeights` or `heliosDiscRasterWeights`
##' @param x starting values for estimated locations x.
##' @param iters number of samples to draw.
##' @return Returns an n x 2 x r array of r samples for the n
##'   locations.
##' @importFrom stats runif
##' @export
heliosDiscSampler <- function(weights,x=NULL,iters=100L) {

  grid <- weights$grid
  model <- weights$model
  W <- weights$W
  m <- nrow(W)
  n <- ncol(W)

  ## Extract model components
  logp.s <- model$logp.s
  logp.p <- model$logp.p
  fixedx <- model$fixedx
  x0 <- model$x0

  if(!is.null(x)) {
    x1 <- x
    q1 <- rep(1,nrow(x))
  } else {
    ## Draw an initial sample
    is <- vapply(seq_len(n),function(k) sample(m,1,prob=W[,k]),0L)
    x1 <- grid[is,]
    q1 <- W[cbind(is,seq_len(n))]
  }
  x1[fixedx] <- x0[fixedx]
  q1[fixedx] <- 1

  ## Allocate chain
  ch <- array(0,c(n,2,iters))

  ## Contribution to logp from each track segment, padded with
  ## fictitious segments at each end
  logp.s1 <- c(0,logp.s(x1),0)

  ## Contributions to the log posterior from the locations
  logp.p1 <- logp.p(x1)


  for(k in seq_len(iters)) {

    ## Propose all x at once, keeping fixed points
    is <- vapply(seq_len(n),function(k) sample(m,1,prob=W[,k]),0L)
    xp <- grid[is,]
    qp <- W[cbind(is,seq_len(n))]
    xp[fixedx] <- x1[fixedx]
    qp[fixedx] <- 1

    for(rb in 1L:2L) {

      ## Modify every second location
      is <- seq.int(rb,n,by=2L)
      x2 <- x1
      q2 <- q1
      x2[is,] <- xp[is,]
      q2[is] <- qp[is]

      ## Contributions to the log posterior from each track
      ## segment, padded with a fictious segment at each end.
      logp.s2 <- c(0,logp.s(x2),0)

      ## Contributions to the log posterior from the locations
      logp.p2 <- logp.p(x2)

      ## The contribution to the log posterior from each location are
      ## the sums of the contributions from the surrounding segments.
      logp1 <- logp.s1[is]+logp.s1[is+1L]+logp.p1[is]
      logp2 <- logp.s2[is]+logp.s2[is+1L]+logp.p2[is]
      logqr <- log(q2[is]) - log(q1[is])


      ## MH rule - compute indices of the accepted points.
      accept <- is[logp2-logp1-logqr > log(runif(length(is)))]
      x1[accept,] <- x2[accept,]
      q1[accept] <- q2[accept]
      logp.p1[accept] <- logp.p2[accept]
      logp.s1[accept] <- logp.s2[accept]
      logp.s1[accept+1L] <- logp.s2[accept+1L]

    }
    ## Store the current state
    ch[,,k] <- x1
  }
  ch
}




##' Compute an approximate track on a discrete spatial grid by
##' maximizing the Helios posterior regularizing by distance
##' travelled.
##'
##' The `heliosDiscMinimalTrack` function seeks the shortest track
##' defined on a discrete spatial grid that is consistent with the
##' observed light, by maximizing the contribution to the posterior
##' from the light data regularized by total distance travelled.
##'
##' This function uses dynamic programming to maximise the log Helios
##' posterior penalized by the total distance travelled, when the
##' track points are constrained to a discrete spatial grid.  This
##' function requires a weights object calculated with
##' [`heliosDiscGridWeights()`] or [`heliosDiscRasterWeights()`] with
##' `log=TRUE`.
##'
##' The `gamma` argument scales the distance penalty used to
##' regularize the posterior.  Larger values of `gamma` result in
##' shorter tracks, at the expense of the accuracy with which the
##' fitted track reproduces the observed light.  For efficiency, the
##' distance penalty is calculated for only the `ngbrs` nearest
##' neighbours of each grid point.
##'
##' @title Minimal distance paths on a discrete grid.
##' @param weights a weights object generated with
##'   `heliosDiscGridWeights` or `heliosDiscRasterWeights`
##' @param gamma regularization parameter controlling the distance cost.
##' @param ngbrs the number of neighbours of each grid point to consider.
##' @return track locations as a two column matrix.
##' @export
heliosDiscMinimalTrack <- function(weights,gamma,ngbrs=50) {

  ## Extract weights and grid
  grid <- weights$grid
  W <- weights$W
  forbid <- weights$model$forbid
  n <- nrow(grid)

  ## Fixed points
  x0 <- weights$model$x0
  fixedx <- weights$model$fixedx

  ## Cache cos(lat) and sin(lat) on the grid
  rad <- pi/180
  coslat <- cos(rad*grid[,2L])
  sinlat <- sin(rad*grid[,2L])

  ## For each grid point calculate distances and indices of near
  ## neighbours
  D <- matrix(0.0,n,ngbrs)
  N <- matrix(0L,n,ngbrs)
  for(i in seq_len(n)) {
    d <- 6378.137*acos(pmin.int(
      coslat*coslat[i]*cos(rad*(grid[,1L]-grid[i,1L]))+sinlat*sinlat[i],1))
    ks <- order(d)[1L:ngbrs]
    N[i,] <- ks
    D[i,] <- d[ks]
  }

  ## Map fixed points to grid
  fixed <- integer(ncol(W))
  for(k in seq_along(fixed)) {
    if(fixedx[k])
      fixed[k] <- which.min(6378.137*acos(pmin.int(
        coslat*cos(rad*x0[k,2L])*
          cos(rad*(grid[,1L]-x0[k,1L]))+sinlat*sin(rad*x0[k,2L]),1)))
  }
  
  ## Path matrix
  K <- matrix(0L,nrow(W),ncol(W))

  ## Tabulate log posterior along paths
  if(fixed[1L]>0)
    maxL <- ifelse(1L:n==fixed[1L],0,forbid)
  else
    maxL <- log(W[,1L])
  

  for(k in 2L:ncol(W)) {

    ## Neither point fixed
    if(fixed[k-1L]==0 && fixed[k]==0) {
      ## Log posterior at each point from each neighbour
      L <- maxL[N]-gamma*D
      ## Select the neighbours that maximise L
      ks <- cbind(1L:n,max.col(L))
      ## Update maximum      
      maxL <- L[ks]+log(W[,k])
      ## Record path
      K[,k-1L] <- N[ks]

    } else {

      ## Both points fixed
      if(fixed[k-1L]>0 && fixed[k]>0) {
        ## Reset maximum
        maxL <- ifelse(1L:n==fixed[k],0,forbid)
        ## Record path
        K[,k-1L] <- fixed[k-1L]
      }
      
      ## Previous point fixed
      if(fixed[k-1L]>0 && fixed[k]==0) {
        ## All distances to fixed point
        kf <- fixed[k-1L]
        d <- 6378.137*acos(pmin.int(
          coslat*coslat[kf]*
            cos(rad*(grid[,1L]-grid[kf,1L]))+sinlat*sinlat[kf],1))
        ## Update maximum
        maxL <- maxL[kf]-gamma*d+log(W[,k])
        ## Record path
        K[,k-1L] <- kf
      }

      ## Current point fixed
      if(fixed[k-1L]==0 && fixed[k]>0) {
        ## All distances to fixed point
        kf <- fixed[k]
        d <- 6378.137*acos(pmin.int(
          coslat*coslat[kf]*
            cos(rad*(grid[,1L]-grid[kf,1L]))+sinlat*sinlat[kf],1))
        ## Record path
        K[,k-1L] <- which.max(maxL-gamma*d)
        ## Reset maximum
        maxL <- ifelse(1L:n==fixed[k],0,forbid)
      }
    }

  }

  ## Backtrace paths
  p <- integer(ncol(W))
  p[ncol(W)] <- k <- which.max(maxL)
  if(!is.finite(maxL[k])) warning("Minimal track inconsistent with data")
  for(k in (ncol(W)-1L):1L)
    p[k] <- K[p[k+1L],k]

  ## Return coords of path
  grid[p,]
}


##' Number of locations
##'
##' A convenience function to determine the number of locations a
##' chain, or set of initial locations or a location summary
##' describe. Assumes `s` is either an array or a list of arrays in
##' which the first dimension corresponds to location, and returns the
##' length of the first dimension.
##'
##' @title Number of locations
##' @param s an array or a list of arrays.
##' @return size of the first dimension of the array.
##' @export
nlocation <- function(s) {
  dim(if(is.list(s)) s[[1L]] else s)[1L]
}


##' Summarize a set of location samples
##'
##' These functions compute various summaries of a sample or list of
##' samples generated by `heliosMetropolis`.
##'
##' These functions accept either a sample from a single mcmc run, or
##' a list of samples from parallel mcmc runs.  When `collapse` is
##' true, multiple samples are merged and single result is returned,
##' otherwise a result is returned for each sample.
##'
##' @rdname locationSummary
##' @title Summaries of Location Samples
##' @param obj the object, or a single chain or a list of parallel
##'   chains generated by `heliosMetropolis`.
##' @param time the times corresponding to the x locations.
##' @param discard number of initial samples to discard.
##' @param alpha coverage of the credible intervals calculated by
##'   `locationSummary`.
##' @param collapse whether to collapse parallel chains to a single
##'   chain
##' @param chains the set of chains to retain, or `NULL`.
##' @return
##' \item{`locationSummary`}{returns a dataframe or a list of
##' dataframes of summary quantities for each location.}
##' \item{`locationMean`}{returns an array or a list of arrays
##' of the means of the samples for each location.}
##' @importFrom stats sd quantile
##' @export
locationSummary <- function(obj,time=NULL,discard=0,alpha=0.95,
                            collapse=TRUE,chains=NULL) {
  summary <- function(s) {
    stat <- function(x) c(mean=mean(x),
                          sd=sd(x),
                          quantile(x,prob=c(0.5,(1-alpha)/2,1-(1-alpha)/2)))
    lon <- t(apply(s[,1L,],1L,stat))
    colnames(lon) <- paste("Lon",colnames(lon),sep=".")
    lat <- t(apply(s[,2L,],1L,stat))
    colnames(lat) <- paste("Lat",colnames(lat),sep=".")
    d <- as.data.frame(cbind(lon,lat))
    if(!is.null(time)) {
      ## Add timing information
      n <- nrow(d)
      if(length(time)==n)
        d <- cbind(Time=time,d)
      else
        d <- cbind(Time1=time[1L:n],Time2=time[2L:(n+1L)],d)
    }
    d
   }

  if(is.list(obj) && !is.null(obj$x)) {
    s <- obj$x
    time <- obj$model$time
  } else
    s <- obj
  
  s <- chainCollapse(s,collapse=collapse,discard=discard,chains=chains)
  if(is.list(s)) lapply(s,summary) else summary(s)
}

##' @rdname locationSummary
##' @export
locationMean <- function(obj,discard=0,collapse=TRUE,chains=NULL) {
  locmean <- function(s) apply(s[,1L:2L,],1L:2L,mean)

  s <- if(is.list(obj) && !is.null(obj$x)) obj$x else obj
  s <- chainCollapse(s,collapse=collapse,discard=discard,chains=chains)
  if(is.list(s)) lapply(s,locmean) else locmean(s)
}


##' Bin locations to form a 2D density image
##'
##' Bins the samples for a sequence of locations to produce 2D array
##' suitable for plotting with `image`.  Optionally, a vector of
##' weights can be provided to differentially weight samples by
##' location.
##'
##' This function accepts either a sample from a single mcmc run, or a
##' list of samples from parallel mcmc runs.  If `collapse` is true,
##' multiple samples are merged and single image is returned,
##' otherwise an image is returned for each sample.
##'
##' @title Location Density Image
##' @param s a single chain or a list of parallel chains generated by
##' `heliosMetropolis`.
##' @param xlim range of the first coordinate.
##' @param ylim range of the second coordinate.
##' @param nx number of cells in the first coordinate.
##' @param ny number of cells in the second coordinate.
##' @param weight weights for each location.
##' @param discard number of initial samples to discard.
##' @param collapse whether to collapse parallel chains to a single chain
##' @param chains the set of chains to retain, or `NULL`.
##' @return A list with elesments
##' \item{`x`}{the x-ordinates that bound the bins}
##' \item{`y`}{the y-ordinates that bound the bins}
##' \item{`W`}{the weighted image.}
##' @export
locationImage <- function(s,xlim,ylim,nx,ny,weight=rep_len(1,dim(s)[1L]),
                           discard=0,collapse=TRUE,chains=NULL) {
  nx <- round(nx)
  ny <- round(ny)
  xbin <- seq.int(xlim[1L],xlim[2L],length.out=nx+1L)
  ybin <- seq.int(ylim[1L],ylim[2L],length.out=ny+1L)

  bin <- function(s) {
    W <- 0
    for(k in 1L:dim(s)[1L]) {
      W <- W+weight[k]*table(
        factor(.bincode(s[k,1L,],xbin),levels=1L:nx),
        factor(.bincode(s[k,2L,],ybin),levels=1L:ny))
    }
    W[W==0] <- NA
    list(x=xbin,y=ybin,W=W)
  }

  s <- chainCollapse(s,collapse=collapse,discard=discard,chains=chains)
  if(is.list(s)) lapply(s,bin) else bin(s)
}




##' Manipulate the samples generated by the Metropolis samplers.
##'
##' These functions provide some basic operations on the samples
##' generated by the Metropolis samplers for Helios.
##'
##' @rdname chainSummary
##' @title Manipulate MCMC samples
##' @param s a single chain or a list of parallel chains generated by
##'   `heliosMetropolis`.
##' @param discard number of initial samples to discard.
##' @param thin rate at which to thin the sample.
##' @param collapse whether to collapse parallel chains to a single
##'   chain
##' @param chains the set of chains to retain, or `NULL`.
##' @return
##' * `chainSummary` returns a summary of the sample
##' * `chainTail` discards the initial samples from each chain
##' * `chainLast` returns the last sample for each location in each chain
##' * `chainCollapse` collapses multiple chains into a single sample
##' * `chainCov` returns the covariance of the parameters location by location
##'   as a pxpxn array.
##' * `chainBcov` returns the joint covariance of the parameters
##'   as an (np)x(np) array.
##' * `chainAcceptance` returns the acceptance rate in the (thinned) chain
##' @export
chainSummary <- function(s) {
  dm <- dim(s[[1]])
  cat("Sample of",
      dm[3L],"from",
      length(s),"chains of",
      dm[2L],"parameters for",
      dm[1L],"locations\n")
}

##' @rdname chainSummary
##' @export
chainTail <- function(s,discard=0,thin=1) {
  tail <- function(s)
    s[,,seq.int(from=1+max(discard,0),to=dim(s)[3L],by=thin)]
  if(!is.list(s)) tail(s) else lapply(s,tail)
}

##' @rdname chainSummary
##' @export
chainLast <- function(s) {
  last <- function(s) s[,,dim(s)[3L]]
  if(!is.list(s)) last(s) else lapply(s,last)
}


##' @rdname chainSummary
##' @export
chainCollapse <- function(s,collapse=TRUE,discard=0,thin=1,chains=NULL) {
  subset <- function(s)
    s[,,seq.int(from=1+max(discard,0),to=dim(s)[3L],by=thin)]
  if(!is.list(s)) {
    if(thin>1 || discard>0)
      s <- subset(s)
  } else {
    if(!is.null(chains)) s <- s[chains]
    if(thin>1 || discard>0) s <- lapply(s,subset)
    if(collapse) {
      dm <- dim(s[[1]])
      s <- array(unlist(s),c(dm[1L:2L],length(s)*dm[3L]))
    }
  }
  s
}



##' @rdname chainSummary
##' @importFrom stats var
##' @export
chainCov <- function(s,discard=0,chains=NULL) {
  s <- chainCollapse(s,collapse=FALSE,discard=discard,chains=chains)

  if(!is.list(s)) {
    V <- apply(s,1L,function(x) var(t(x)))
  } else {
    dm <- dim(s[[1]])
    V <- apply(array(
      unlist(lapply(s, function(s) apply(s,1L,function(y) var(t(y))))),
      c(dm[c(2L,2L,1L)],length(s))),
      1L:3L,mean)
  }
  V
}



##' @rdname chainSummary
##' @importFrom stats var
##' @export
chainBcov <- function(s,discard=0,chains=NULL) {
  bcov <- function(s) {
    dm <- dim(s)
    dim(s) <- c(prod(dm[1L:2L]),dm[3L])
    var(t(s))
  }

  s <- chainCollapse(s,collapse=FALSE,discard=discard,chains=chains)
  if(is.list(s))
    apply(simplify2array(lapply(s,bcov)),1L:2L,mean)
  else
    bcov(s)
}

##' @rdname chainSummary
##' @export
chainAcceptance <- function(s,collapse=FALSE,chains=NULL) {
  rate <- function(s)
    mean(apply(s,1,function(x) mean(rowMeans(x[,-1L]-x[,-ncol(x)]!=0))))

  s <- chainCollapse(s,collapse=FALSE,chains=chains)
  r <- if(is.list(s)) lapply(s,rate) else rate(s)
  if(collapse & is.list(r)) do.call(mean,r) else r
}



##' Convert to Coda objects
##'
##' Convert samples generated by `helioMetropolis` or to a \pkg{coda}
##' object.
##'
##' @title Export to Coda
##' @param s a list of chains generated by `heliosMetropolis`
##' @return a \pkg{coda} object.
##' @importFrom coda mcmc mcmc.list
##' @export
chainCoda <- function(s) {
  coda <- function(s) {
      dm <- dim(s)
      dim(s) <- c(prod(dm[1L:2L]),dm[3L])
      nms <- c("Lon","Lat")
      if(dm[2L]>2)
          nms <- c("Lon","Lat",paste0("P",seq.int(length.out=dm[2L]-2)))
      nms <- as.vector(t(outer(nms,1L:dm[1L],paste,sep=".")))
      rownames(s) <- nms
      mcmc(t(s))
  }
  if(is.list(s)) {
      if(length(s)==1) coda(s[[1L]]) else do.call(mcmc.list,lapply(s,coda))
  } else {
      coda(s)
  }
}




##' Construct a sampler to draw multivariate Normal deviates.
##'
##' Construct a sampler that draws Multivariate Normal deviates with
##' covariances determined by `S` and mean determined by its
##' first argument.
##'
##' The `mvnorm` function constructs a function that generates `n`
##' independent Multivariate Normal deviates, where the covariance of
##' each deviate is specified by `S` which must be either a `m`x`m`
##' covariance matrix or an `n`x`m`x`m` array of covariance matrices.
##'
##' The `bmvnorm` function constructs a function that generates `n`
##' correlated Multivariate Normal deviates, where the joint
##' covariance is specified by `S` which must be a `nm`x\`nm`
##' covariance matrix (as generated by `chainBcov`).
##'
##' @title Multivariate Normal Samples
##' @param S a covariance matrix or an array of covariance matrices.
##' @param s a scale factor applied to S.
##' @param n number of deviates to draw.
##' @param m dimension of each deviate.
##' @param tol minimum allowable variance.
##' @return A function that draws bivariate Normal deviates with mean
##'   given by its first argument.
##' @importFrom stats rnorm
##' @export
mvnorm <- function(S,s=1,n=1,tol=1.0E-6) {

  ## Fault tolerant cholesky
  fchol <- function(V) {
    diag(V) <- pmax.int(diag(V),tol)
    tryCatch(chol(V),
             error=function(e) {
               d <- diag(V)
               V <- V/4
               diag(V) <- pmax.int(d/2,tol)
               chol(V)
             })
  }

  m <- dim(S)[1L]
  if(length(dim(S))==2) {
    S <- array(s*fchol(S),c(m,m,n))
  } else {
    n <- dim(S)[3L]
    for(k in 1L:n) {
      S[,,k] <- s*fchol(S[,,k])
    }
  }
  S <- aperm(S,c(1L,3L,2L))
  dim(S) <- c(m*n,m)

  function(mu) {
    z <- rnorm(m*n)*S
    dim(z) <- c(m,m*n)
    z <- colSums(z)
    dim(z) <- c(n, m)
    mu+z
  }
}

##' @rdname mvnorm
##' @importFrom stats rnorm
##' @export
bmvnorm <- function(S,m,s=1) {
  S <- chol(s*S)
  n <- nrow(S)/m

  function(mu) {
    z <- rnorm(m*n)%*%S
    dim(z) <- c(n,m)
    mu+z
  }
}

