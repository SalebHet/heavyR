#' mvnpdf
#'
#' Function compute all the values of the density evaluated at the n points in a single call of the function
#'
#' @export
#' @import mvtnorm
#' @param x Matrix of size NxP
#' @param mean a vector of means
#' @param varcovM a variance-covariance matrix
#' @param Log a logical parameter, with default value to \code{TRUE}
#'
#' @return \code{mat} matrix x
#' @return \code{vect}  a vector of length n of the multivariate normal distribution density values at those points.
#'
#' @examples
#' M1<-matrix(runif(36),nrow=6)
#' res <- heavyR::mvnpdf(M1)
#' res$mat
#' res$vect

mvnpdf <- function(x,mean = rep(0,nrow(x)),varcovM = diag(nrow(x)),Log = TRUE){
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))
  const <- p/2 * log(2*pi) - 0.5 * LogDetvarcovM

  res <- NULL
  for (j in 1:ncol(x)) { #idée opti lapply()
    yj <- - const -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    res <- c(res, yj)
  }

  if(!Log){
    res <- exp(res)
  }
  return(list(mat = x,vect = res))
}
