# SRS.R

#' @name SRS
#'
#' @title Simple random sampling.
#'
#' @description This function invokes base::sample() to draw a random sample using
#' a user specified random seed.
#'
#' @details This function was written by Phil Davies.
#'
#' @param seed The random seed to be used to draw the current sample.
#' @param total_rows The total number of rows in our input dataset.
#' @param sample_size The number of rows wanted in our random sample.
#'
#' @return A random sample.

#' @examples
#' # Create a random sample with a seed of 99 ----------------------------------
#' spbal::SRS(seed = 99, total_rows = 100, sample_size = 20)
#'
#' # Create a random sample with a seed of 42 ----------------------------------
#' spbal::SRS(seed = 42, total_rows = 100, sample_size = 20)
#'
#' # Create a random sample with a seed of 99 ----------------------------------
#' spbal::SRS(seed = 99, total_rows = 100, sample_size = 25)
#'
#' @export
SRS <- function(seed = 511, total_rows = 0, sample_size = 0) {

  # validate our parameters.
  validate_parameters("seed", base::c(seed))
  validate_parameters("total_rows", base::c(total_rows))
  validate_parameters("sample_size", base::c(sample_size))

  # ensure sample_size < total_rows
  if(sample_size >= total_rows){
    base::stop("spbal(SRS) Parameter sample_size must be less than total_rows.")
  }

  samp_pts <- NULL
  base::set.seed(seed)
  samp_pts <- base::sample(x      = total_rows,
                          size    = sample_size,
                          replace = FALSE,
                          prob    = NULL)
  return (samp_pts)
}

