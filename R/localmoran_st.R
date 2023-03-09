# Based on spdep::localmoran by Roger Bivand Roger.Bivand@nhh.no
# Copyright 2001-18 by Roger Bivand, 2021 Jeff Sauer and Levi Wolf (conditional code), 2023 Nils Tjaden (spatio-temporal adaption)
# package version: 1.2.2
# license: GPL-2 | GPL-3 [expanded from: GPL (≥ 2)]
#

#' @importFrom stats na.fail naresid
#'
localmoran_st <- function(x, listw, zero.policy=NULL, na.action=na.fail,
                          conditional=TRUE, alternative = "two.sided",
                          mlvar=TRUE, spChk=NULL, adjust.x=FALSE) {
  ## input checks
  #stopifnot(is.vector(x))
  if (!inherits(listw, "listw")) # check if the weight list is a weights list
    stop(paste(deparse(substitute(listw)), "is not a listw object"))
  if (is.null(zero.policy))      # if no zero policy is defined, get it from the spdep options
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
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
  # Ii
  #res[,1] <- (z/s2) * lz
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
