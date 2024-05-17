# getSample.R

#' @name getSample
#'
#' @title Extract a sample of a specified size from a master sample.
#'
#' @description A description of this useful function.
#'
#' @author Phil Davies.
#'
#' @param shapefile A MULTIPOINT or POINT object from where to take the sample.
#' @param n The number of sample points to return.
#' @param randomStart Whether a spatially balanced sample will be randomly drawn from
#' the frame or not. Default value is FALSE.
#' @param strata to be added
#' @param stratum The name of a column in the dataframe attached to shapefile that defines
#' the strata of interest.
#'
#' @return A list containing the following variable:
#'
#' \itemize{
#' \item \code{sample} The sample from the shapefile POINTS.
#' }
#'
#' @examples
#' # Draw a spatially balanced sample of n = 25 from a Halton Frame over Gates --
#'
#' # Use the North Carolina shapefile supplied in the sf R package.
#' shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
#' shp_gates <- shp_file[shp_file$NAME == "Gates",]
#'
#' # Vertically aligned master sample bounding box.
#' bb <- spbal::BoundingBox(shapefile = shp_gates)
#'
#' set.seed(511)
#' result7 <- spbal::HaltonFrame(shapefile = shp_gates,
#'                               J = c(6, 4),
#'                               boundingbox = bb)
#' Frame <- result7$hf.pts.shp
#'
#' # Get the first 25 sites from a B = (2^6) * (3^4) Halton Frame (62,208 grid
#' # points covering Gates).
#' n_samples <- 25
#' FrameSample <-getSample(shapefile = Frame,
#'                         n = n_samples)
#' FrameSample <- FrameSample$sample
#' FrameSample
#'

#' @export
getSample <- function(shapefile,
                      n,
                      randomStart = FALSE,
                      strata = NULL,
                      stratum = NULL){

  # need to validate what function created shapefile. Only accept output from the spbal::HaltonFrame
  # function and spbal::BAS.
  if(!base::attr(shapefile, "spbal") %in% base::c("HaltonFrame", "BAS")){
    base::stop("spbal(getSample) unsupported function for this spbal shapefile.")
  }

  # strata and stratum must be both NULL or both not NULL.
  if((base::is.null(strata) && !base::is.null(stratum))|(!base::is.null(strata) && base::is.null(stratum))){
    base::stop("spbal(getSample) strata and stratum must both have valid values or both be NULL.")
  }
  # are we dealing with a stratified sample? default is not.
  stratification <- FALSE
  if(!base::is.null(strata) && !base::is.null(stratum)){
    stratification <- TRUE
    # random start not suported if stratification.
    randomStart = FALSE
  } else {
    # only check n when not using stratification.
    # can accept input of either a MULTIPOINT object or a POINT object.
    # either $sample or $hf.pts.shp from the HaltonFrame function.
    validate_parameters("n", base::c(n))
  }

  # Get the geometry type
  geometry_type <- sf::st_geometry_type(shapefile)

  if(!base::any(geometry_type == "POINT")){
    if(!base::any(geometry_type == "MULTIPOINT")){
      base::stop("spbal(getSample) Supplied shapefile must contain POINT or MULTIPOINT geometries.")
    }
  }

  # Check if the geometry type is POINT or MULTIPOINT
  is_point <- base::all(geometry_type == "POINT")
  is_multipoint <- base::all(geometry_type == "MULTIPOINT")

  # if a MULTIPOINT object need to make a POINT object first. should be from $hf.pts.shp.
  if(is_multipoint){
    # get data for the Halton frame.
    #base::message("spbal(getSample) is_multipoint.")
    hf_pts <- sf::st_cast(shapefile, "POINT")
    hf_pts <- sf::st_as_sf(hf_pts)
    hf_pts$ID <- base::seq(1, base::length(hf_pts$x))
    hf_pts$spbalSeqID <- seq(1, base::length(hf_pts$x))
    shapefile <- hf_pts
    # now we can make our sample based on randomStart.
  }

  # need to ensure n is less than the number of points in shapefile.

  # if is_point TRUE then shapefile from $sample
  if(randomStart){
    # points will already be sorted and each have a valid $spbalSeqID
    #base::message("spbal(getSample) is_randomStart.")
    duplicated_pts <- base::rbind(shapefile, shapefile)
    random_start_point <- base::sample(1:base::length(duplicated_pts$ID), 1)
    sample_indices <- base::seq(random_start_point, (random_start_point + n) - 1, 1)
    sample <- duplicated_pts[sample_indices,]
  } else {
    # points will already be sorted and each have a valid $spbalSeqID
    # lets sort on $spbalSeqID to be sure...
    shp_sorted <- shapefile[base::order(shapefile$spbalSeqID), ]
    if(stratification){
      # then we will already have provided n samples.
      stratified_sample <- shp_sorted[shp_sorted[[stratum]] %in% strata,]
      #if(n > base::length(stratified_sample$spbalSeqID)){
      #  msg <- "spbal(getSample) Warning - n (%s) exceeds number of points in strata (%s). Returning %s points."
      #  msgs <- base::sprintf(msg, n, strata, base::length(stratified_sample$spbalSeqID))
      #  base::message(msgs)
      n <- base::length(stratified_sample$spbalSeqID)
      #}
      # get our n sample points.
      sample <- stratified_sample[1:n,]
    } else {
      # not a random start so just return the first n sample points.
      #base::message("spbal(getSample) is_not_randomStart and not stratified.")
      if(n > base::length(shp_sorted$spbalSeqID)){
        msg <- "spbal(getSample) Warning - n (%s) exceeds number of points in shapefile. Returning %s points."
        msgs <- base::sprintf(msg, n, base::length(shp_sorted$spbalSeqID))
        base::message(msgs)
        n <- base::length(shp_sorted$spbalSeqID)
      }
      # get our n sample points.
      sample <- shp_sorted[1:n,]
    }
  }

  # assign the spbal attribute to the sample being returned, i.e. the function that created it.
  base::attr(sample, "spbal") <- "getSample"

  # package up object to be returned.
  result <- base::list(sample = sample)
  return(result)
}

