# findFirstStudyRegionPoint.R

#' @name findFirstStudyRegionPoint
#'
#' @title Get a randomly chosen Halton point from within the study area and the associated seeds.
#'
#' @description This function repeatedly calls function spbal::getBASSample
#' to generate the Halton frame sample. This function selects the first point at random from those
#' points in the study area. This point and the seeds used to generate the sample are returned to
#' the caller.
#'
#' @author This function was written by Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) of the study area(s).
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL and will
#' be defined randomly using function \code{generateUVector}.
#' @param bb Bounding box which defines the Master Sample. A bounding box must be
#' supplied.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing three variables:
#'
#' \itemize{
#' \item \code{seeds} The u1 and u2 seeds used to generate the first point.
#' \item \code{k} The index of the first point in the initial sample.
#' }
#'
# 1. Set J1 = 4 and J2 = 3.
# 2. Generate B = 2^J1 x 3^J2 points from a random-start Halton sequence H
#    with a random seed (u1, u2).
# 3. Find points from H in the study area. Call this set S. If S is empty,
#    increment J1 and J2 and go to step 2.
# 4. Randomly choose a point from S. Let xk be this point where k is the 'site index'
#    (I think that's what we call it).
# 5. Set the seed to (u1 + k - 1, u2 + k - 1).
# 6. Re-number ID by subtracting k (re-generate sample using seeds for 5 - first point must also be ID=1)

# For example, let (u1, u2) = (1, 5) and S = {x2, x6, x7}.
# If x6 is randomly chosen, then the new seed is (1 + 6 - 1, 5 + 6 - 1) = (6, 10)
# (the sixth point in H).

# The only difference is that the random-start Halton sequence must be length B.
#' @keywords internal
findFirstStudyRegionPoint <- function(shapefile, bb, seeds, verbose = FALSE){
  # must not be called without seeds! (also checked in getBASSample).
  if(base::is.null(seeds)){
    msg <- "spbal(findFirstStudyRegionPoint) The seeds parameter must not be NULL."
    msgs <- base::sprintf(msg)
    base::stop(msgs)
  }

  # Initialise variables.
  J <- base::c(4, 3)
  bases <- base::c(2, 3)
  crs <- sf::st_crs(shapefile)

  # default number of sample points to find.
  n <- (bases[1]^J[1]) * (bases[2]^J[2])

  pts_in_intersection <- 0
  call.getBASSample.cnt <- 0

  while(pts_in_intersection < 1){
    # shapefile, bb, n, seeds
    call.getBASSample.cnt <- call.getBASSample.cnt + 1
    result <- getBASSample(shapefile = shapefile, bb = bb, n = n, seeds = seeds)
    diff_ <- result$sample
    seeds <- result$seed

    # find number of points within our study area/bb intersection.
    pts_in_intersection <- base::length(diff_$SiteID)
    n <- n * 2
  }

  if(verbose){
    msg <- "spbal(findFirstStudyRegionPoint) Needed %s call(s) to locate first study area point."
    msgs <- base::sprintf(msg, call.getBASSample.cnt)
    base::message(msgs)
  }

  # select a point in the study area at random.
  base::set.seed(seeds[1] + seeds[2])
  k <- base::sample(pts_in_intersection, 1)
  k <- diff_$SiteID[k]

  # select our random first point. # return SiteID = 1, know what k is.
  #first.pt <- diff_ #[1,]

  if(verbose){
    msg <- "spbal(findFirstStudyRegionPoint) Random point selected: %s."
    msgs <- base::sprintf(msg, k)
    base::message(msgs)
  }

  result <- base::list(seeds = seeds,
                       k     = k)
  return(result)
}

