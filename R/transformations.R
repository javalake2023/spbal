# transformations.R

#' @name rot
#'
#' @title Generate a rotation matrix for rotating objects later.
#'
#' @description Generate a rotation matrix for rotating objects later.
#'
#' @author This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param a radians of rotation.
#'
#' @return Matrix
#'
#' @keywords internal
rot <- function(a){
  base::matrix(base::c(base::cos(a), base::sin(a), -base::sin(a), base::cos(a)), 2, 2)
}


#' @name rotate.scale.coords
#'
#' @title Scale and rotate points from the unit square to a defined projection.
#'
#' @description Given some coordinates on \[0,1)x\[0,1), shift and scale them to the bounding box, and then rotate
#' them given the bounding box rotation defined by the Master Sample.
#'
#' @author This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample.
#'
#' @param coords Output from RSHalton() to be converted to the spatial surface of interest.
#' @param bb Special shape file defining the bounding box with attributes for centroid and rotation.
#' @param back Boolean for whether or not the rotation is back to the original rotated bounding box.
#'
#' @return sf spatial points with projection defined in bb.
#'
#' @keywords internal
rotate.scale.coords <- function(coords, bb, back = TRUE){

  coords <- coords[, 2:3]
  theta <- ifelse(back, -1, 1) * base::attr(bb, "rotation")	# Rotate backwards
  cntrd <- base::attr(bb, "centroid")

  bb.bounds <- sf::st_bbox(bb)
  bb.scale <- base::diag(2) * (bb.bounds[3:4] - bb.bounds[1:2])

  coords.tmp <- sf::st_multipoint(coords, dim = "XY")
  coords <- sf::st_geometry(coords.tmp)
  coords.scale <- coords * bb.scale + bb.bounds[1:2]
  coords.rot <- (coords.scale - cntrd) * rot(theta) + cntrd

  sf::st_crs(coords.rot) <- sf::st_crs(bb)
  return(coords.rot)
}

