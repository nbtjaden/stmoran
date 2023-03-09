#' Calculate CORT
#'
#' Calculates first order temporal correlation coefficient (CORT) between two
#' time series. This is based on equation 3 in Gao et al. (2019), which is in
#' turn from page 10 of Chouakria & Nagabhushan (2007)
#'
#' @author Nils Tjaden <n.b.tjaden@@utwente.nl>
#' @references Gao, Y. et al. (2019) ‘Measuring spatio-temporal autocorrelation
#'   in time series data of collective human mobility’, Geo-spatial Information
#'   Science, 22(3), pp. 166–173. doi:10.1080/10095020.2019.1643609.
#' @param X a vector, containing the data values of the first time series
#' @param Y a vector, containing the data values of the second time series
#'
#' @return an atomic vector, containing the CORT value [-1,1]
#'
#' @examples
#' data("EuStockMarkets")
#' CORT(X=EuStockMarkets[ ,1], Y=EuStockMarkets[ ,2])
#'
CORT <- function(X, Y){
  n <- length(X)
  CORT <- sum((X[2:n] - X[1:(n-1)]) * (Y[2:n] - Y[1:(n-1)])) / (sqrt(sum((X[2:n] - X[1:(n-1)])^2)) * sqrt(sum((Y[2:n] - Y[1:(n-1)])^2)))
  if(is.nan(CORT)){
    warning("NaN produced. This most likely means that X, Y, or both have no variation (i.e. all values are identical), which leads to a division by 0. Sorry, CORT() can not handle this properly at the moment.")
  }
  return(CORT)
}

#' Tuning function (phi) for CORT
#'
#' Tunes the first order correlation coefficient from [-1,1] to (0,2).
#' This is based on equation 5 in Gao et al. (2019).
#'
#' @author Nils Tjaden <n.b.tjaden@@utwente.nl>
#' @references Gao, Y. et al. (2019) ‘Measuring spatio-temporal autocorrelation
#'   in time series data of collective human mobility’, Geo-spatial Information
#'   Science, 22(3), pp. 166–173. doi:10.1080/10095020.2019.1643609.
#' @param x a vector, containing the result of CORT()
#'
#' @return a vector, containing the tuned value
#'
#' @examples
#' data("EuStockMarkets")
#' tune_CORT(CORT(X=EuStockMarkets[ ,1], Y=EuStockMarkets[ ,2]))
#'
#' plot(tune_CORT(seq(from=-1, to=1, by=.01))~seq(from=-1, to=1, by=.01), type="l", xlab="CORT value", ylab="tuned value")
#'
tune_CORT <- function(x){
  2 / (1 + exp(2*x))
}

#' Calculate magnitude deviation
#'
#' The wording in Gao et al. (2019) describing is slightly obscure. They make no
#' direct reference to any of the methods used in Chouakria & Nagabhushan
#' (2007), and indeed it appears to be a modification of the Euclidian distance
#' method that does not enforce positive values. The last paragraph of section
#' 2.3 of Gao et al. (2019) explains that: "The quantity deviation determines
#' whether, and how much, the autocorrelation is positive or negative; then the
#' similarity of the temporal variance adjusts the degree of the correlation."
#'
#' @author Nils Tjaden <n.b.tjaden@@utwente.nl>
#' @references Gao, Y. et al. (2019) ‘Measuring spatio-temporal autocorrelation
#'   in time series data of collective human mobility’, Geo-spatial Information
#'   Science, 22(3), pp. 166–173. doi:10.1080/10095020.2019.1643609.
#' @param X a vector, containing the data values of the first time series
#' @param Y a vector, containing the data values of the second time series
#'
#' @return an atomic vector, containing the magnitude deviation
#'
#' @examples
#' data("EuStockMarkets")
#' magnitude_deviation(X=EuStockMarkets[ ,1], Y=EuStockMarkets[ ,2])
#'
magnitude_deviation <- function(X, Y){
  sum(X) - sum(Y) # note: sum(X - Y) is equivalent but slower
}

#' Calculate time series deviation
#'
#' Measure for the deviation between two time series. This can be used to
#' estimate the difference between a specific time series and the mean time
#' series.
#'
#' @author Nils Tjaden <n.b.tjaden@@utwente.nl>
#' @references Gao, Y. et al. (2019) ‘Measuring spatio-temporal autocorrelation
#'   in time series data of collective human mobility’, Geo-spatial Information
#'   Science, 22(3), pp. 166–173. doi:10.1080/10095020.2019.1643609.
#' @param X a vector, containing the data values of the first time series
#' @param Y a vector, containing the data values of the second time series
#'
#' @return an atomic vector, containing the deviation between the two time
#'   series
#'
#' @examples
#' data("EuStockMarkets")
#' ts_deviation(X=EuStockMarkets[ ,1], Y=EuStockMarkets[ ,2])
#'
ts_deviation <- function(X, Y){
  #tune_CORT(CORT(X, Y)) * (sum(X) - sum(Y))
  tune_CORT(CORT(X, Y)) * magnitude_deviation(X, Y)
}
