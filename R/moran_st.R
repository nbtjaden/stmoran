# Based on spdep::moran by Roger Bivand Roger.Bivand@nhh.no
# Copyright 2001-18 by Roger Bivand, 2023 Nils Tjaden (spatio-temporal adaption)
# package version: 1.2.2
# license: GPL-2 | GPL-3 [expanded from: GPL (≥ 2)]

#' Calculate spatio-temporal Moran's I
#'
#' Calculates a spatio-temporal version of (global) Moran's I as proposed by Gao
#' et al. (2019). Called by `moran.test_st`.
#'
#' @export
#' @author Roger Bivand <Roger.Bivand@@nhh.no> (original code of spdep::moran)
#' @author Nils Tjaden <n.b.tjaden@@utwente.nl> (adjustments for spatio-temporal
#'   use)
#' @references Gao, Y. et al. (2019) ‘Measuring spatio-temporal autocorrelation
#'   in time series data of collective human mobility’, Geo-spatial Information
#'   Science, 22(3), pp. 166–173. doi:10.1080/10095020.2019.1643609.
#'
#' @param x a numeric matrix (or data.frame) of observations. Columns represent
#'   space, rows represent time, so that each column contains a time series
#'   measured at one location. The number of columns must equal the length of
#'   the neighbors list in `listw`.
#' @param listw a `listw` object (see documentation for `spdep`)
#' @param n number of zones
#' @param S0 global sum of weights
#' @param zero.policy default NULL, use global option value; if TRUE assign zero
#'   to the lagged value of zones without neighbours, if FALSE assign NA
#' @param NAOK 	if 'TRUE' then any 'NA' or 'NaN' or 'Inf' values in x are passed
#'   on to the foreign function. If 'FALSE', the presence of 'NA' or 'NaN' or
#'   'Inf' values is regarded as an error.
#' @returns  A list: `I` spatio-temporal Moran's I `K` sample kurtosis of x
#' @examples
moran_st <- function (x, listw, n, S0, zero.policy = NULL, NAOK = FALSE){
  if (is.null(zero.policy))
    zero.policy <- spdep::get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  w <- spdep::listw2mat(listw = listw)
  if (n1 != length(x))
    stop("objects of different length")
  xx <- rowMeans(x, na.rm=NAOK)
  z <- rep(x = NA, times=n1)
  for(i in 1:n1){
    z[i] <- ts_deviation(X=x[,i], Y=xx)
  }

  zz <- sum(z^2, na.rm = NAOK)
  K <- (length(x) * sum(z^4, na.rm = NAOK))/(zz^2)
  lz <- spdep::lag.listw(listw, z, zero.policy = zero.policy, NAOK = NAOK)

  I <- (n/S0) * ((sum(z * apply(X=w, MARGIN = 1, FUN=function(x){sum(x*z)})))/zz)

  res <- list(I = I, K = K)
  res
}

#' Moran's I test for spatial autocorrelation (spatio-temporal)
#'
#' This is a 1:1 copy of `spdep::moran.test` with minor adjustments for use with
#' `moran_st`. See `?spdep::moran.test` for more information.
#' @export
#'
#' @author Roger Bivand <Roger.Bivand@@nhh.no> (original code of
#'   spdep::moran.test)
#' @author Nils Tjaden <n.b.tjaden@@utwente.nl> (minor adjustments for
#'   spatio-temporal use)
#' @references Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion
#' @param x a numeric matrix (or data.frame) of observations. Columns represent
#'   space, rows represent time, so that each column contains a time series
#'   measured at one location. The number of columns must equal the length of
#'   the neighbors list in `listw`.
#' @param listw a `listw` object (see documentation for `spdep`)
#' @param randomisation variance of I calculated under the assumption of
#'   randomisation, if FALSE normality
#' @param zero.policy default NULL, use global option value; if TRUE assign zero
#'   to the lagged value of zones without neighbours, if FALSE assign NA
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of greater (default), less or two.sided.
#' @param rank logical value - default FALSE for continuous variables, if TRUE,
#'   uses the adaptation of Moran's I for ranks suggested by Cliff and Ord
#'   (1981, p. 46)
#' @param na.action a function (default `na.fail`), can also be na.omit or
#'   na.exclude - in these cases the weights list will be subsetted to remove
#'   NAs in the data. It may be necessary to set zero.policy to TRUE because
#'   this subsetting may create no-neighbour observations. Note that only
#'   weights lists created without using the glist argument to `nb2listw` may be
#'   subsetted. If `na.pass` is used, zero is substituted for NA values in
#'   calculating the spatial lag
#' @param adjust.n default TRUE, if FALSE the number of observations is not
#' adjusted for no-neighbour observations, if TRUE, the number of observations
#' is adjusted
#' @param spChk	should the data vector names be checked against the spatial
#'   objects for identity integrity, TRUE, or FALSE, default NULL to use
#'   `get.spChkOption()`
#' @param drop.EI2	default FALSE, if TRUE, emulate CrimeStat <= 4.02
#' @returns A list with class htest containing the following components:
#' @returns `statistic` the value of the standard deviate of Moran's I.
#' @returns `p.value` the p-value of the test.
#' @returns `estimate`	the value of the observed Moran's I, its expectation and
#'   variance under the method assumption.
#' @returns `alternative` a character string describing the alternative
#'   hypothesis.
#' @returns `method`	a character string giving the assumption used for
#'   calculating the standard deviate.
#' @returns `data.name` a character string giving the name(s) of the data.
#' @examples
moran.test_st <- function(x, listw, randomisation=TRUE, zero.policy=NULL,
                          alternative="greater", rank = FALSE, na.action=na.fail, spChk=NULL,
                          adjust.n=TRUE, drop.EI2=FALSE) {
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
                                            "is not a listw object"))
  if (!is.numeric(unlist(x))) stop(paste(deparse(substitute(x)),
                                         "is not a numeric vector"))
  if (is.null(zero.policy))
    zero.policy <- spdep::get.ZeroPolicyOption()
  stopifnot(is.logical(zero.policy))
  if (is.null(spChk)) spChk <- spdep::get.spChkOption()
  if (spChk && !spdep::chkIDs(t(x), listw))
    stop("Check of data and weights ID integrity failed")
  xname <- deparse(substitute(x))
  wname <- deparse(substitute(listw))
  NAOK <- deparse(substitute(na.action)) == "na.pass"
  x <- na.action(x)
  na.act <- attr(x, "na.action")
  if (!is.null(na.act)) {
    subset <- !(1:length(listw$neighbours) %in% na.act)
    listw <- subset(listw, subset, zero.policy=zero.policy)
  }
  n <- length(listw$neighbours)
  if (n != length(x)) stop("objects of different length")

  wc <- spdep::spweights.constants(listw, zero.policy=zero.policy,
                            adjust.n=adjust.n)
  S02 <- wc$S0*wc$S0
  res <- moran_st(x, listw, wc$n, wc$S0, zero.policy=zero.policy,
                  NAOK=NAOK)
  I <- res$I
  K <- res$K
  if (rank) K <- (3*(3*wc$n^2 -7))/(5*(wc$n^2 - 1))
  EI <- (-1) / wc$n1
  if(randomisation) {
    VI <- wc$n*(wc$S1*(wc$nn - 3*wc$n + 3) - wc$n*wc$S2 + 3*S02)
    tmp <- K*(wc$S1*(wc$nn - wc$n) - 2*wc$n*wc$S2 + 6*S02)
    if (tmp > VI) warning("Kurtosis overflow,\ndistribution of variable does not meet test assumptions")
    VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
    if (!drop.EI2) VI <- (VI - EI^2)
    if (VI < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  } else {
    VI <- (wc$nn*wc$S1 - wc$n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
    if (!drop.EI2) VI <- (VI - EI^2)
    if (VI < 0) warning("Negative variance,\ndistribution of variable does not meet test assumptions")
  }
  ZI <- (I - EI) / sqrt(VI)
  statistic <- ZI
  names(statistic) <- "Moran I statistic standard deviate"
  if (alternative == "two.sided"){
    PrI <- 2 * stats::pnorm(abs(ZI), lower.tail=FALSE)
  } else {
    if (alternative == "greater"){
      PrI <- stats::pnorm(ZI, lower.tail=FALSE)
    } else {PrI <- stats::pnorm(ZI)}
  }
  if (!is.finite(PrI) || PrI < 0 || PrI > 1) {
    warning("Out-of-range p-value: reconsider test arguments")
  }
  vec <- c(I, EI, VI)
  names(vec) <- c("spatio-temporal Moran I statistic", "Expectation", "Variance")
  method <- paste("spatio-temporal Moran I test under", ifelse(randomisation,
                                                               "randomisation", "normality"))
  data.name <- paste(xname, ifelse(rank,
                                   "using rank correction",""), "\nweights:",
                     wname, ifelse(is.null(na.act), "", paste("\nomitted:",
                                                              paste(na.act, collapse=", "))),
                     ifelse(adjust.n && isTRUE(any(sum(spdep::card(listw$neighbours) == 0L))),
                            "n reduced by no-neighbour observations\n", ""),
                     ifelse(drop.EI2, "EI^2 term dropped in VI", ""), "\n")
  res <- list(statistic=statistic, p.value=PrI, estimate=vec,
              alternative=alternative, method=method, data.name=data.name)
  if (!is.null(na.act)) attr(res, "na.action") <- na.act
  class(res) <- "htest"
  res
}
