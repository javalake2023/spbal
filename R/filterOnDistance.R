# filterOnDistance.R

#' @name filterOnDistance
#'
#' @title Filter sample using a minimum distance.
#'
#' @description The input parameter minRadius >= 0 is the minimum distance between any two
#' points in the sample. My idea is to apply this condition to the points in the over-sample,
#' result$overSample. Let's call these points x1, x2, ..., xB. Create a new set S = (x1).
#' Starting from x1, we check if dist(S,x2) > minRadius. If it is, add x2 to S. For x3, we check
#' if dist(S,x3) > minRadius, where dist is the smallest distance from a point in S to x3
#' (single linkage distance). If dist(S,x3) > minRadius, add x3 to S. Continue until you reach xB.
#'
#' The distances are calculated as great circles over an oblate spheroid and the units are meters.
#'
#' @author Phil Davies.
#'
#' @details Key points:
#'
#' \itemize{
#'   \item \code{result$minRadius} is nonempty (it always contains x1). Hence, if the user chooses a crazy minRadius, they get one point.
#'   \item \code{result$minRadius} is a subset of result$overSample.
#'   \item The number of points in result$minRadius is random. That's fine!
#'   \item If they want n points and result$minRadius has less than n, too bad! They can reduce minRadius and/or increase the iterations parameter.
#'   \item If they want a sample with the minimum radius property, they use:
#'    \itemize{
#'             \item \code{smp <- result$minRadius}
#'             \item \code{sample <- smp[1:n,]}
#'    }
#' }
#'
#' @param overSample A HIP sample.
#' @param minRadius The minimum distance between any two points in the sample.
#'
#' @return S The set of points that are more than minRadius from each other.
#'
#' @keywords internal

filterOnDistance <- function(overSample, minRadius){

  # initialise set S.
  S <- NULL
  # if a valid minRadius value has been specified then calculate set S.
  # probably need to place in a function at some stage.
  if(!base::is.null(minRadius)){
    # Initialize the set S with the first point
    S <- overSample[1, , drop = FALSE]
    # code to calculate distance and check against minRadius
    for (i in 2:base::length(overSample$spbalSeqID)) {
      # Calculate distance from current point to all points in S
      distances <- sf::st_distance(overSample[i, , drop = FALSE], S)

      # If the minimum distance is greater than minRadius, add to S
      if("units" %in% class(distances)){
        if (base::min(distances) > units::set_units(minRadius, "m")) {
          S <- base::rbind(S, overSample[i, , drop = FALSE])
        }
      } else {
        if (base::min(distances) > minRadius) {
          S <- base::rbind(S, overSample[i, , drop = FALSE])
        }
      }
    }
  } else {
    S <- overSample
  }
  return(S)
}

