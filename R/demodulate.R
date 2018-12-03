#' demodulate
#'
#' Function to apply complex demodulation, a mathematical time series analysis method, to seasonal infectious disease data.
#'
#' The method of complex demodulation is in three parts. First, the data is multiplied by a demodulating function. Second, a filter is applied (here, a centred moving average). Third, the amplitude and phase are extracted from the filtered, demodulated data. The output is a list of two objects: the amplitude, and the phase.
#'
#' @param x Data. It needs to be in vector format.
#' @param p Dominant periodicity of the data.
#'
#' @return Returns a list of three objects: the input data, the amplitude, and the phase.
#' @references
#' Hogan, A.B., Glass, K. and Anderssen, R.S. (2017) Complex demodulatuion: a novel #' time series method for analysing infectious diseases. \emph{The ANZIAM Journal} \strong{59}(1) 51-60.
#'
#'
#' @seealso \code{\link{synthdata}}
#'
#' @examples
#' x <- synthdata(seq(1,120), 12, 0.2, 0.2)
#' demodulate(x, 12)
#'
demodulate <- function(x,p) {
  if (is.vector(x) == F){
    stop('x needs to be a vector')
  }
  n <- length(x)  # dimensions of data
  N <- n[1]       # data length
  tvec <- seq(1,n[1],1) # time vector
  q <- p/2
  f0 <- 1/p # frequency
  i <- sqrt(as.complex(-1)) # define imaginary i
  s <- exp((-2)*pi*i*f0*tvec) # demodulator
  loopstart <- (q+1)
  loopend <- (N-q)
  X <- as.vector(x*s)
  Y <- vector(mode="complex", length=N)

  is.even <- function(x) x %% 2 == 0
  if (is.even(length(x)) == F) {
    for (j in c(loopstart:loopend)){
      Y[j] <- (1 / p) * (0.5 * X[(j - q)] + 0.5 * X[(j + q)] + sum(X[(j - q + 1):(j + q - 1)]))
    }} else {
      for (j in c(loopstart:loopend)) {
        Y[j] <- (1 / p) * (sum(X[(j - q):(j + q-1)]))
      }
    }

  amp <- 2*Mod(Y) # Complex modulus
  phase <- Im(log(Y))/(2*pi)

  return(list(data = x, amp = amp, phase = phase))
}

#-----------------------------------------------
