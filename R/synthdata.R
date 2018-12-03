#' synthdata
#'
#'Generates synthetic biennial epidemic data.
#'
#' @param t Number of time points for which to generate data
#' @param p Period of the dominant annual cycle
#' @param a Degree of alterating peak sizes, in the range [0,1]
#' @param b Degree of delay in alternating peaks, in the range [0,1]
#'
#' @return The output will be a synthetic timeseries of length t.
#' @references
#' Hogan, A.B., Glass, K. and Anderssen, R.S. (2017) Complex demodulatuion: a novel #' time series method for analysing infectious diseases. \emph{The ANZIAM Journal} \strong{59}(1) 51-60.
#'
#' @seealso \code{\link{demodulate}}
#'
#' @examples
#' synthdata(seq(1,120), 12, 0.3, 0.2)
#'

synthdata <- function(t, p, a = 0, b = 0){
  t <- 1:t
  x <- (1 + a*cos(pi*t/p))*(cos(pi*t/p + b*cos(pi*t/p)))^2 + abs(rnorm(length(t), mean = 1, sd=0.1))
  return(x)
}
