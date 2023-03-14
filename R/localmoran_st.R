# Based on spdep::localmoran by Roger Bivand Roger.Bivand@nhh.no
# Copyright 2001-18 by Roger Bivand, 2021 Jeff Sauer and Levi Wolf (conditional code), 2023 Nils Tjaden (spatio-temporal adaption)
# package version: 1.2.2
# license: GPL-2 | GPL-3 [expanded from: GPL (≥ 2)]
#

#' Calculate spatio-temporal Moran's I
#'
#' Calculates a spatio-temporal version of local Moran's I (a.k.a
#' spatio-temporal LISA) as proposed by Gao et al. (2019).
#'
#' @importFrom stats na.fail naresid
#' @export
#' @author Roger Bivand <Roger.Bivand@@nhh.no> (original code of spdep::moran)
#' @author Jeff Sauer (conditional code)
#' @author Levi Wolf (conditional code)
#' @author Nils Tjaden <n.b.tjaden@@utwente.nl> (adjustments for spatio-temporal
#'   use)
#' @references Gao, Y. et al. (2019) ‘Measuring spatio-temporal autocorrelation
#'   in time series data of collective human mobility’, Geo-spatial Information
#'   Science, 22(3), pp. 166–173. doi:10.1080/10095020.2019.1643609.
#' @references Anselin, L. 1995. Local indicators of spatial association,
#'   Geographical Analysis, 27, 93–115
#' @references Getis, A. and Ord, J. K. 1996 Local spatial statistics: an
#'   overview. In P. Longley and M. Batty (eds) Spatial analysis: modelling in a
#'   GIS environment (Cambridge: Geoinformation International), 261–277
#' @references Sokal, R. R, Oden, N. L. and Thomson, B. A. 1998. Local Spatial
#'   Autocorrelation in a Biological Model. Geographical Analysis, 30. 331–354
#' @references Bivand RS, Wong DWS 2018 Comparing implementations of global and
#'   local indicators of spatial association. TEST, 27(3), 716–748
#'   doi:10.1007/s11749-018-0599-x
#' @references Sauer, J., Oshan, T. M., Rey, S., & Wolf, L. J. 2021. The
#'   Importance of Null Hypotheses: Understanding Differences in Local Moran’s
#'   under Heteroskedasticity. Geographical Analysis. doi:10.1111/gean.12304
#' @references Bivand, R. (2022), R Packages for Analyzing Spatial Data: A
#'   Comparative Case Study with Areal Data. Geographical Analysis, 54(3),
#'   488-518. doi:10.1111/gean.12319
#'
#' @param x a numeric matrix (or data.frame) of observations. Columns represent
#'   space, rows represent time, so that each column contains a time series
#'   measured at one location. The number of columns must equal the length of
#'   the neighbors list in `listw`.
#' @param listw a `listw` object (see documentation for `spdep`)
#' @param zero.policy default NULL, use global option value; if TRUE assign zero
#'   to the lagged value of zones without neighbours, if FALSE assign NA
#' @param na.action a function (default `na.fail`), can also be na.omit or
#'   na.exclude - in these cases the weights list will be subsetted to remove
#'   NAs in the data. It may be necessary to set zero.policy to TRUE because
#'   this subsetting may create no-neighbour observations. Note that only
#'   weights lists created without using the glist argument to `nb2listw` may be
#'   subsetted. If `na.pass` is used, zero is substituted for NA values in
#'   calculating the spatial lag
#' @param conditional default TRUE: expectation and variance are calculated
#'   using the conditional randomization null (Sokal 1998 Eqs. A7 & A8).
#'   Elaboration of these changes available in Sauer et al. (2021). If FALSE:
#'   expectation and variance are calculated using the total randomization null
#'   (Sokal 1998 Eqs. A3 & A4).
#' @param conditional	default TRUE: expectation and variance are calculated
#'   using the conditional randomization null (Sokal 1998 Eqs. A7 & A8).
#'   Elaboration of these changes available in Sauer et al. (2021). If FALSE:
#'   expectation and variance are calculated using the total randomization null
#'   (Sokal 1998 Eqs. A3 & A4).
#' @param alternative	a character string specifying the alternative hypothesis,
#'   must be one of greater, less or two.sided (default).
#' @param mlvar	default TRUE: values of local Moran's I are reported using the
#'   variance of the variable of interest (sum of squared deviances over n), but
#'   can be reported as the sample variance, dividing by (n-1) instead; both are
#'   used in other implementations.
#' @param spChk	should the data vector names be checked against the spatial
#'   objects for identity integrity, TRUE, or FALSE, default NULL to use
#'   get.spChkOption()
#' @param adjust.x default FALSE, if TRUE, x values of observations with no
#'   neighbours are omitted in the mean of x

#' @returns `Ii_st` local moran statistic
#' @returns `E.Ii_st` expectation of local moran statistic
#' @returns `Var.Ii_st` variance of local moran statistic
#' @returns `Z.Ii_st` standard deviate of local moran statistic
#' @returns `Pr()_st` p-value of local moran statistic using pnorm()
#' @returns In addition, an attribute data frame "quadr" indicating the
#'   High-High etc. quadrants
#'
localmoran_st <- function(x, listw, zero.policy=NULL, na.action=na.fail,
                          conditional=TRUE, alternative = "two.sided",
                          mlvar=TRUE, spChk=NULL, adjust.x=FALSE) {
  ## input checks
  if (!inherits(listw, "listw")) # check if the weight list is a weights list
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(zero.policy))      # if no zero policy is defined, get it from the spdep options
    zero.policy <- spdep::get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  alternative <- match.arg(alternative, c("two.sided", "greater", "less")) # verify that alternative hypotheses is specified correctly
  if (!is.null(attr(listw$neighbours, "self.included")) &&
      attr(listw$neighbours, "self.included"))
    stop("Self included among neighbours")
  if (is.null(spChk)) spChk <- spdep::get.spChkOption() # check IDs for data integrity?
  if (spChk && !spdep::chkIDs(t(x), listw))
    stop("Check of data and weights ID integrity failed")
  if (!is.numeric(unlist(x)))
    stop(paste(deparse(substitute(x)), "is not a numeric vector"))

  # handling NAs
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action") # will be null if nothing was done to existing NA's or no NA's exist

  rn <- attr(listw, "region.id")
  if (!is.null(na.act)) { # if anything was done to the data when checking NA's, adjust weights list accordingly
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy=zero.policy)
    excl <- class(na.act) == "exclude"
  }

  n <- length(listw$neighbours)
  if (n != length(x))stop("Different numbers of observations")
  res <- matrix(nrow=n, ncol=5)

  # alternative hypothesis
  if (alternative == "two.sided"){
    Prname <- "Pr(z != E(Ii))"
  } else {
    if (alternative == "greater") {
      Prname <- "Pr(z > E(Ii))"
    } else {
      Prname <- "Pr(z < E(Ii))"
    }
  }

  colnames(res) <- paste0(c("Ii", "E.Ii", "Var.Ii", "Z.Ii", Prname), "_st")

  # global mean time series
  if (adjust.x) { # omit values of observations with no neighbours in the mean of x
    nc <- spdep::card(listw$neighbours) > 0L
    xx <- rowMeans(x[nc], na.rm=NAOK)
  } else {
    xx <- rowMeans(x, na.rm=NAOK)
  }

  # calculate time series deviations
  # (deviation from mean)
  #z <- x - xx
  z <- rep(x = NA, times=n)
  for(i in 1:n){
    z[i] <- ts_deviation(X=x[,i], Y=xx)
  }

  # determining quadrants.
  lz <- spdep::lag.listw(listw, z, zero.policy=zero.policy, NAOK=NAOK)
  lbs <- c("Low", "High")
  quadr <- interaction(cut(z, c(-Inf, 0, Inf), labels=lbs),
                       cut(lz, c(-Inf, 0, Inf), labels=lbs), sep="-")

  # calculate variance of data (deviation from mean)
  if (mlvar) { # divide by n or n-1?
    if (adjust.x) { # omit x values of observations with no neighbors?
      s2 <- sum(z[nc]^2, na.rm=NAOK)/sum(nc)
    } else {
      s2 <- sum(z^2, na.rm=NAOK)/n
    }
  } else {
    if (adjust.x) {
      s2 <- sum(z[nc]^2, na.rm=NAOK)/(sum(nc)-1)
    } else {
      s2 <- sum(z^2, na.rm=NAOK)/(n-1)
    }
  }

  # calculate spatio-temporal local Moran's I
  Ii <- rep(x = NA, times=n)
  for(i in 1:n){
    neighbours <- listw$neighbours[[i]]
    weights <- listw$weights[[i]]

    Ii[i] <- (n*z[i] * sum(weights * z[neighbours])) / sum(z^2)
  }
  res[,1] <- Ii

  # calculate expectation (E.i)
  Wi <- sapply(listw$weights, sum) # sum up the weights of all neighbours per location, results in either 0 or 1. See Sokal 1998 section 3.1
  if (conditional){	# use conditional randomization null rather than total randomization null (Sokal 1998 eq. A7 vs. A3). See also Sauer et al. (2021).
    m2 <- sum(z * z) / n
    res[, 2] <- -(z ** 2 * Wi) / ((n - 1) * m2)
  } else {
    res[, 2] <- -Wi / (n-1)
  }

  # calculate variance (Var.Ii)
  if (mlvar)  {# divide by n or n-1?
    if (adjust.x) {# omit x values of observations with no neighbors?
      b2 <- (sum(z[nc]^4, na.rm=NAOK)/sum(nc))/(s2^2)
    } else {
      b2 <- (sum(z^4, na.rm=NAOK)/n)/(s2^2)
    }
  } else {
    if (adjust.x) {
      b2 <- (sum(z[nc]^4, na.rm=NAOK)/(sum(nc)-1))/(s2^2)
    } else {
      b2 <- (sum(z^4, na.rm=NAOK)/(n-1))/(s2^2)
    }
  }

  Wi2 <- sapply(listw$weights, function(x) sum(x^2)) #See Sokal 1998 section 3.1
  A <- (n-b2) / (n-1)
  B <- (2*b2 - n) / ((n-1)*(n-2))
  if (conditional){ # use conditional randomization null rather than total randomization null (Sokal 1998 eq. A8 vs. A4). See also Sauer et al. (2021).
    res[, 3] <- ((z / m2) ** 2 *
                   (n / (n - 2)) *
                   (Wi2 - (Wi ** 2 / (n - 1))) *
                   (m2 - (z ** 2 / (n - 1))))
  } else { # conditional=FALSE
    res[,3] <- A*Wi2 + B*(Wi^2 - Wi2) - res[,2]^2
  }


  # calculate Z.Ii (standard deviate of local moran statistic)
  res[,4] <- (res[,1] - res[,2]) / sqrt(res[,3])

  # calculate p-values
  if (alternative == "two.sided") {
    pv <- 2 * stats::pnorm(abs(res[,4]), lower.tail=FALSE)
  } else if (alternative == "greater") {
    pv <- stats::pnorm(res[,4], lower.tail=FALSE)
  } else {
    pv <- stats::pnorm(res[,4])
  }
  res[,5] <- pv

  # adjust residuals for missing values
  if (!is.null(na.act) && excl) {
    res <- naresid(na.act, res)
  }
  if (!is.null(rn)) rownames(res) <- rn
  attr(res, "call") <- match.call()
  if (!is.null(na.act)) attr(res, "na.action") <- na.act
  class(res) <- c("localmoran_st", class(res))
  attr(res, "quadr") <- data.frame(quadr=quadr)

  res
}
