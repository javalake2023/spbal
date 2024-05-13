# BoundingBox.R

#' @name BoundingBox
#'
#' @title Create a bounding box for a study region.
#'
#' @description Randomly generate a seed from 10,000 possible values in right now 2 dimensions.
#' Note that in van Dam-Bates et al. (2018) we required that the random seed falls into main
#' object shape, such as one of the islands in New Zealand, or within marine environment for
#' BC west coast. However, with a random rotation, we are able to ignore that detail. If this
#' function is used without a random rotation, we recommend running it until
#' the first master sample point does indeed fall within the largest scale of the master sample use.
#'
#' @author This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later ported to this package, spbal.
#'
#' @param shapefile Spatial feature that defines the boundary of the area to define a bounding
#' box over.
#' @param d Dimension of the new Master Sample, at this stage we only work with d=2.
#' @param rotate Boolean of whether or not to randomly rotate the bounding box. This parameter
#' is not supported at this time.
#' @param verbose Print the rotation and random seed when it is generated.
#'
#' @return bounding box for a study area.
#'
#' @examples
#' # Create a bounding box for the Gates, North Carolina study area -------------
#' # Use the North Carolina shapefile supplied in the sf R package.
#' shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
#' shp_gates <- shp_file[shp_file$NAME == "Gates",]

#' # Vertically aligned master sample bounding box.
#' bb <- spbal::BoundingBox(shapefile = shp_gates)
#' bb
#'
#' @export
BoundingBox <- function(shapefile, d = 2, rotate = FALSE, verbose = FALSE){

  # validate the shapefile parameter.
  bb_geometry <- sf::st_geometry_type(shapefile)
  if (!base::all(bb_geometry %in% base::c("MULTIPOLYGON", "POINT", "POLYGON", "MULTIPOINT"))){
    msg <- "spbal(BoundingBox) Unsupported geometry in shapefile, %s."
    msgs <- base::sprintf(msg, bb_geometry)
    base::stop(msgs)
  }

  # Just work with sf objects.
  if (base::class(shapefile)[1] != "sf"){
    msg <- "spbal(BoundingBox) Shapefile does not have class of sf, %s, converting to sf object."
    msgs <- base::sprintf(msg, base::class(shapefile)[1])
    base::message(msgs)
    shp <- sf::st_as_sf(shapefile)
  }

  # We always use base (2, 3) in this version of spbal.
  base <- base::c(2, 3, 5)[1:d]

  # rotate is not supported at this time, will always set theta to 0.
  if(rotate) {
    theta <- stats::runif(1, -base::pi, base::pi)
  }else{
    theta <- 0
  }

  # create bounding box.
  bb <- sf::st_as_sfc(sf::st_bbox(shapefile))
  cntrd <- sf::st_centroid(bb)
  base::attr(bb, "rotation") = theta
  base::attr(bb, "centroid") = sf::st_coordinates(cntrd)
  sf::st_crs(bb) <- sf::st_crs(shapefile)

  # assign the spbal attribute to the object being returned, i.e. the function that created it.
  base::attr(bb, "spbal") <- "BoundingBox"

  return(bb)
}

