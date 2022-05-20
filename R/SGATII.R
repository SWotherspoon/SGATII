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
  cumsum(c(wrapLon(lon[1],lmin),wrapLon(diff(lon))))



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
##' @title Solar Zenith Angle
##' @param sun list of solar time and declination computed by
##' `solar`.
##' @param lon vector of longitudes.
##' @param lat vector latitudes.
##' @return A vector of solar zenith angles (degrees) for the given
##' locations and times.
##' @seealso [solar()]
##' @examples
##' ## Approx location of Sydney Harbour Bridge
##' lon <- 151.211
##' lat <- -33.852
##' ## Solar zenith angle for noon on the first of May 2000
##' ## at the Sydney Harbour Bridge
##' s <- solar(as.POSIXct("2000-05-01 12:00:00","EST"))
##' zenith(s,lon,lat)
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
##' refraction.
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
                            1735+elev*(-518.2+elev*(103.4+elev*(-12.79+elev*0.711))),-20.772/te)))
  ## Corrected Zenith [90-AG]
  zenith-r/3600
}




##' @rdname refracted
##' @importFrom stats uniroot
##' @export
unrefracted <- function(zenith)
  uniroot(function(x) refracted(x)-zenith,c(zenith,zenith+2))



##' Estimate time of sunrise or sunset for a given location given the
##' approximate solar time of twilight
##'
##' Solar declinat`ion and equation of time vary slowly over the day,
##' and so the values of the Solar declination and equation of time at
##' sunrise/sunset can be caclulated approximately if an approximate
##' time of sunrise/sunset is known. The sun's hour angle and hence
##' sunrise/sunset for the required zenith can then be calculated from
##' these approximations.
##'
##' Note this function returns the time of twilight in solar time.
##' @title Solar Time of Sunrise and Sunset
##' @param solar output of `solar` for approximate times of
##' twilight.
##' @param lon vector of longitudes.
##' @param lat vector of latitudes.
##' @param rise logical vector indicating whether to compute rise or
##' set.
##' @param zenith the solar zenith angle that defines twilight.
##' @return a vector of twilight times in solar time (degrees)
##' @seealso [twilight()]
##' @export
twilightSolartime <- function(solar,lon,lat,rise,zenith=96) {
  rad <- pi/180
  cosz <- cos(rad*zenith)
  cosHA <- (cosz-sin(rad*lat)*solar$sinSolarDec)/(cos(rad*lat)*solar$cosSolarDec)
  ## Compute the sun's hour angle from its declination for this location
  hourAngle <- ifelse(rise,360,0)+ifelse(rise,-1,1)*suppressWarnings(acos(cosHA)/rad)
  ## Solar time of sunrise at this zenith angle, lon and lat
  #(hourAngle+180-lon)%%360
  #360*(solar$solarTime%/%360)+solarTime
  solarTime <- (hourAngle+180-lon)%%360
  (solarTime-solar$solarTime+180)%%360-180+solar$solarTime
}


##' Estimate time of sunrise or sunset for a given day and location
##'
##' `twilight` uses an iterative algorithm to estimate times of
##' sunrise and sunset.
##'
##' Note that these functions return the twilight that occurs on the
##' same date GMT as `tm`, and so sunset may occur before sunrise,
##' depending upon latitude.
##'
##' Solar declination and equation of time vary slowly over the day,
##' and so the values of the Solar declination and equation of time at
##' sunrise/sunset are well approximated by their values at 6AM/6PM
##' local time. The sun's hour angle and hence sunrise/sunset for the
##' required zenith can then be caclulates from these approximations.
##' The calculation is then repeated using the approximate
##' sunrise/sunset times to derive more accurate values of the Solar
##' declination and equation of time and hence better approximations
##' of sunrise/sunset.  The process is repreated and is accurate to
##' less than 2 seconds within 2 or 3 iterations.
##'
##' It is possible that sunrise or sunset does occur for a given date
##' and location. When `closest` is `FALSE`, the twilight returned on
##' or before the (UTC) date of `tm`.  When `closest` is `TRUE`,
##' `twilight` attempts to return the twilight closest to the input
##' time `tm`.
##'
##' `sunrise` and `sunset` are simple wrappers for `twilight`.
##' @title Times of Sunrise and Sunset
##' @param tm vector of approximate times of twilight.
##' @param lon vector of longitudes.
##' @param lat vector of latitudes.
##' @param rise logical vector indicating whether to compute rise or
##'   set.
##' @param zenith the solar zenith angle that defines twilight.
##' @param iters number of iteratve refinements made to the initial
##'   approximation.
##' @param closest if `TRUE`, attempt to find the twilight closest to
##'   `tm`.
##' @return a vector of twilight times.
##' @examples
##' ## Approx location of Santa Barbara
##' lon <- -119.7022
##' lat <- 34.4191
##' ## Sunrise and sunset for 8th April 2013 at Santa Barbara
##' day <- as.POSIXct("2013-04-08","GMT")
##' sunrise(day,lon,lat)
##' sunset(day,lon,lat)
##' @export
twilight <- function(tm,lon,lat,rise,zenith=96,iters=3,closest=FALSE) {

  ## Compute date
  date <- as.POSIXlt(tm)
  date$hour <- date$min <- date$sec <- 0
  date <- as.POSIXct(date,"GMT")

  lon <- (lon+180)%%360-180
  ## GMT equivalent of 6am or 6pm local time
  twl <- date+240*(ifelse(rise,90,270)-lon)
  ## Iteratively improve estimate
  for(k in seq_len(iters)) {
    s <- solar(twl)
    s$solarTime <- s$solarTime%%360
    solarTime <- 4*twilightSolartime(s,lon,lat,rise,zenith)-s$eqnTime
    twl <- date+60*solarTime
  }

  if(closest) {
    delta <- (as.numeric(tm)-as.numeric(twl))/3600
    off <- double(length(delta))
    off[delta > 12] <- 86400
    off[delta < -12] <- -86400
    twl <- twilight(tm+off,lon,lat,rise,zenith,iters,FALSE)
  }

  twl
}

##' @rdname twilight
##' @export
sunrise <- function(tm,lon,lat,zenith=96,iters=3,closest=FALSE)
  twilight(tm,lon,lat,rise=TRUE,zenith=zenith,iters=iters,closest=closest)

##' @rdname twilight
##' @export
sunset <- function(tm,lon,lat,zenith=96,iters=3,closest=FALSE)
  twilight(tm,lon,lat,rise=FALSE,zenith=zenith,iters=iters,closest=closest)



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

  6378.137*acos(pmin.int(cosx2[-n]*cosx2[-1L]*cos(rad*(x[-1L,1L]-x[-n,1L]))+sinx2[-n]*sinx2[-1L],1))
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
              cosx2[-n]*sinx2[-1L]-sinx2[-n]*cosx2[-1L]*cos(rad*(x[-1L,1L]-x[-n,1L])))/rad
  ## Difference bs and fold difference into [-180,180)
  wrapLon(bs[1-n]-bs[-1])
}



##' Simulate zenith angles, times and locations of twilight along a
##' specified track.
##'
##' Given times, longitudes and latitudes that specify a template
##' track, `zenithSimulate` interpolates the template onto the new
##' times specified by `tm.out` and computes the solar zenith angle at
##' each point along the new track. Given a dataframe generated by
##' `zenithSimulate`, `twilightSimulate` computes times and locations
##' of sunrise and sunset based on the simulated zenith angles. The
##' `twilightPerturb` adds a given vector of errors (in minutes) to
##' the twilights in a dataframe generated by `twilightSimulate`, in
##' such a way that a positive error causes sunrise to occur later and
##' sunset earlier.
##'
##' @title Solar Zenith and Twilight Simulation
##' @param tm vector of times that specify the template track.
##' @param lon vector of longitude that specify the template track.
##' @param lat vector of latitude that specify the template track.
##' @param tm.out vector of times to which the template is resampled.
##' @param dfz a dataframe generated with `zenithSimulate`.
##' @param zenith the solar zenith angle that defines twilight.
##' @return `zenithSimulate` returns a data frame with
##' components
##' \item{`Date`}{times along the simulated track}
##' \item{`Lon`}{longitudes along the simulated track}
##' \item{`Lat`}{latitudes along the simulated track}
##' \item{`Zenith`}{zenith angles along the simulated track}
##' `twilightSimulate` returns a data frame of twilights with
##' components
##' \item{`Twilight`}{times of twilight}
##' \item{`Rise`}{is this a sunrise}
##' \item{`Lon`}{longitude at twilight}
##' \item{`Lat`}{latitude at twilight}
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



##' @rdname zenithSimulate
##' @export
twilightSimulate <- function(dfz,zenith=96) {

  n <- nrow(dfz)

  ## Compute indexes for sunrise and sunset
  sr.k <- which(dfz$Zenith[-n] >= zenith & dfz$Zenith[-1L] < zenith)
  ss.k <- which(dfz$Zenith[-n] < 96 & dfz$Zenith[-1L] >= 96)
  ## Interleave sunrise and sunset
  ord <- order(c(sr.k,ss.k))
  k <- c(sr.k,ss.k)[ord]
  rise <- rep(c(T,F),c(length(sr.k),length(ss.k)))[ord]
  ## Interpolation weights
  w <- (zenith-dfz$Zenith[k])/(dfz$Zenith[k+1L]-dfz$Zenith[k])

  ## Interplated times and locations of twilight
  data.frame(Twilight=dfz$Date[k] + w*(as.vector(dfz$Date[k+1L])-as.vector(dfz$Date[k])),
             Rise=rise,
             Lon=dfz$Lon[k] + w*(dfz$Lon[k+1L]-dfz$Lon[k]),
             Lat=dfz$Lat[k] + w*(dfz$Lat[k+1L]-dfz$Lat[k]))
}




##' Phaethon Model Structures.
##'
##' Phaethon requires a model structure that describes the model being
##' fitted. This function generate basic model structures that should
##' provide a suitable starting point for most analyses.
##'
##' The `phaethonModel` function constructs a model structure that
##' combines the `phaethonLightModel` with a behavioural model that
##' assume the average speed of travel between successive locations is
##' Gamma distributed.
##'
##' The phaeton light model assumes the light sensor cannot detect
##' light when it is dark, but may fail to detect light when it is
##' light because the sensor is obscured.  The model assumes the
##' length of time the sensor can be obscured is exponentially
##' distributed, with rate `alpha`.
##'
##' The speed of travel is assumed Gamma distributed.  The shape and
##' rate of the Gamma distribution of speeds is specified by the
##' `beta` parameter. If `beta` is a two element vector the same shape
##' and rate are applied across all segments, but if `beta` is a two
##' column matrix the shape and rate can be specified on a segment by
##' segment basis. By default, the speed of travel is calculated based
##' on the time intervals between the twilights (in hours), but the
##' intervals of time actually available for travel can be specified
##' directly with the `dt` argument.
##'
##' The `forbid` argument represents the likelihood of light being
##' observed when it is dark.  When `forbid=-Inf` the initial
##' locations `x0` must be consistent with the observed light - there
##' can be no light observed when it is dark on the implied track.
##' Setting `forbid` to a large negative value relaxes this
##' constraint, making the observation of light in the dark unlikely
##' but not impossible.
##'
##' @title Phaethon Model Structures
##' @param date POSIXct vector of times for day/nights observations
##' @param light vector of day/night observations
##' @param tm times at which to fit locations
##' @param x0 initial estimates of fitted locations.
##' @param alpha rate parameter for sensor obscuration.
##' @param beta parameters of the behavioural model.
##' @param logp.p function to evaluate any additional contribution to
##'   the log posterior from the estiamted locations
##' @param fixedx logical vector indicating which twilight locations
##'   to hold fixed.
##' @param dt (optional) vector of time intervals (hours) for speed
##'   calculation.
##' @param zenith the solar zenith angle that defines twilight.
##' @param forbid the log likelihood for forbidden light observations
##' @importFrom stats dgamma
##' @return a list with components \item{`logp`}{function to evaluate
##'   the contributions to the log posterior from the t}
##'   \item{`logpz`}{function to evaluate the contributions to the log
##'   posterior from the prior for the z locations}
##'   \item{`estelle.logpb`}{function to evaluate contribution to the
##'   log posterior from the behavioural model for estelle.}
##'   \item{`stella.logpb`}{function to evaluate contribution to the
##'   log posterior from the behavioural model for stella.}
##'   \item{`residuals`}{function to evaluate the twilight model
##'   residuals.}  \item{`fixedx`}{a logical vector indicating which
##'   locations should remain fixed.}  \item{`x0`}{an array of initial
##'   twilight locations.}  \item{`time`}{the twilight times.}
##'   \item{`rise`}{the sunrise indicators.}  \item{`group`}{the
##'   grouping vector.}
##' @export
phaethonModel <- function(date,light,tm,x0,
                          alpha,beta,
                          logp.p=function(x) rep.int(0L,nrow(x)),
                          fixedx=FALSE,dt=NULL,zenith=96,forbid=-Inf) {

  ## Times (hours) between observations
  if(is.null(dt))
    dt <- diff(as.numeric(tm)/3600)

  ## Ensure beta is always a matrix
  if(!is.matrix(beta)) beta <- t(beta)


  ## Contribution to log posterior for the light data
  logp.x <- phaethonLightModel(date,light,tm,alpha,forbid,zenith)

  logp.b <- function(x) {
    spd <- pmax.int(trackDist(x), 1e-06)/dt
    dgamma(spd,beta[,1L],beta[,2L],log=TRUE)
  }

  logp <- function(x) logp.p(x)+logp.x(x)+logp.b(x)

  list(
    ## Suggested starting points
    tm=tm,
    x0=x0,
    fixedx=fixedx,
    ## Data
    date=date,
    light=light,
    ## Components of the posterior
    logp.p=logp.p,
    logp.b=logp.b,
    logp.x=logp.x,
    logp=logp)
}

##' @rdname phaethonModel
##' @export
phaethonLightModel <- function(date,light,tm,alpha,forbid,zenith) {

  ## Convert twilights to solar time.
  s <- solar(date)

  segments <- cut(date,tm,include.lowest=TRUE)
  ls <- split(light > 0,segments)

  function(x) {
    ## Interpolate
    lon <- approx(tm,x[,1],date,rule=2)$y
    lat <- approx(tm,x[,2],date,rule=2)$y
    ## Calculate where should be day
    zs <- split(zenith(s,lon,lat) < zenith,segments)
    ## Calculate costs per segment
    mapply(function(z,l) if(any(!z & l)) forbid else alpha*sum(z & !l),zs,ls)
  }
}


##' Metropolis sampler for Phaethon
##'
##' These functions draw samples form posterior for the Phaethon
##' model by the Metropolis algorithm.
##'
##' @title Metropolis Samplers
##' @param model a model structure as generated by `thresholdModel`.
##' @param proposal function for drawing proposals for x.
##' @param x0 Starting values for twilight locations x.
##' @param iters number of samples to draw.
##' @param thin rate at which to thin samples.
##' @param chains number of chains to sample.
##' @param verbose report progress at prompt?
##' @return If there are r samples drawn for each of q chains of p
##'   parameters at n locations, Phaethon will return a list
##'   containing \item{`model`}{the model structure} \item{`x`}{a list
##'   of n x p x r arrays of twilight locations from the q chains}
##' @importFrom stats runif
##' @importFrom utils flush.console
##' @export
phaethonMetropolis <- function(model,proposal,x0=NULL,
                              iters=1000L,thin=10L,chains=1L,
                              verbose=interactive()) {

  ## Initialize x
  if(is.null(x0)) x0 <- model$x0
  ## Expand starting values for multiple chains
  x0 <- rep(if(is.list(x0)) x0 else list(x0),length.out=chains)

  ## Number of locations
  n <- nrow(x0[[1]])
  ## Number of parameters
  m <- ncol(x0[[1]])

  ## Extract model components
  logp <- model$logp
  fixedx <- model$fixedx

  ## Lists of chains
  ch.xs <- vector(mode="list",chains)

  ## PARALLEL - parallelise this loop
  for(k1 in 1:chains) {
    ## Allocate chain
    ch.x <- array(0,c(n,m,iters))

    ## Initialize
    x1 <- x0[[k1]]
    ## Drop dimnames for speed
    dimnames(x1) <- NULL

    ## Contribution to logp from each track segment, padded with
    ## fictious segments at each end
    logp.s1 <- c(0,logp(x1),0)


    k2 <- 0
    if(verbose) {
      cat("iter ",sprintf("%6d",k2))
      flush.console()
    }

    for(k2 in 1:iters) {

      if(verbose && k2%%10==0) {
        cat("\b\b\b\b\b\b");
        cat(sprintf("%6d",k2));
        flush.console()
      }

      for(k3 in 1:thin) {

        ## Propose all x at once, keeping fixed points, and calculate
        ## contribution to the log posterior
        xp <- proposal(x1)
        xp[fixedx,] <- x1[fixedx,]

        for(rb in 1:2) {

          ## Modify every second location
          is <- seq.int(rb,n-1L,by=2L)
          x2 <- x1
          x2[is,] <- xp[is,]

          ## Contributions to the log posterior from each track
          ## segment, padded with a fictious segment at each end.
          logp.s2 <- c(0,logp(x2),0)

          ## The contribution to the log posterior from each location are
          ## the sums of the contributions from the surrounding segments.
          logp1 <- logp.s1[is]+logp.s1[is+1L]
          logp2 <- logp.s2[is]+logp.s2[is+1L]

          ## MH rule - compute indices of the accepted points.
          accept <- is[logp2-logp1 > log(runif(length(is)))]
          x1[accept,] <- x2[accept,]
          logp.s1[accept] <- logp.s2[accept]
          logp.s1[accept+1L] <- logp.s2[accept+1L]
        }
      }
      ## Store the current state
      ch.x[,,k2] <- x1
    }
    ch.xs[[k1]] <- ch.x
    if(verbose) cat("\n")
  }
  list(model=model,x=ch.xs)
}


