# generateUVector.R

#' @name generateUVector
#'
#' @title Generate a vector of two random seeds.
#'
#' @description This function generates two seeds, u1 and u2, in the range 0 to 2^11 and 0 to 3^7
#' respectively. These are returned to the caller in the form of a vector. This is for internal
#' use only.
#'
#' @author Phil Davies.
#'
#' @return A vector containing two seeds, u1 and u2.
#'
#' @keywords internal
generateUVector <- function(){
  u1 <- base::floor(stats::runif(1, 0, 2^11))
  u2 <- base::floor(stats::runif(1, 0, 3^7))
  seeds <- base::c(u1, u2)
  return(seeds)
}

