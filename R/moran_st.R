# Based on spdep::moran by Roger Bivand Roger.Bivand@nhh.no
# Copyright 2001-18 by Roger Bivand, 2023 Nils Tjaden (spatio-temporal adaption)
# package version: 1.2.2
# license: GPL-2 | GPL-3 [expanded from: GPL (≥ 2)]


moran_st <- function (x, listw, n, S0, zero.policy = NULL, NAOK = FALSE){
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
  stopifnot(is.logical(zero.policy))
  n1 <- length(listw$neighbours)
  w <- spdep::listw2mat(listw = listw)
  #x <- c(x)
  if (n1 != length(x))
    stop("objects of different length")
  xx <- rowMeans(x, na.rm=NAOK)
  #z <- x - xx
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


moran.test_st <- function(x, listw, randomisation=TRUE, zero.policy=NULL,
                          alternative="greater", rank = FALSE, na.action=na.fail, spChk=NULL,
                          adjust.n=TRUE, drop.EI2=FALSE) {
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
                                            "is not a listw object"))
  if (!is.numeric(unlist(x))) stop(paste(deparse(substitute(x)),
                                         "is not a numeric vector"))
  if (is.null(zero.policy))
    zero.policy <- get("zeroPolicy", envir = .spdepOptions)
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
