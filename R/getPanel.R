# getPanel.R

#' @name contains_feature
#'
#' @title Check if the sf object contains a specified feature.
#'
#' @description Used to check if a simple file object contains a feature. This
#' is an internal only function.
#'
#' @author Phil Davies.
#'
#' @param sf_object Simple file object that we want to verify if it contains
#' a feature called feature_name.
#' @param feature_name The feature name we want to find in the simple file
#' object sf_object.
#'
#' @return Returns TRUE if the simple file object sf_object contains the feature
#' feature_name. Otherwise FALSE is returned.
#'
#' @keywords internal
contains_feature <- function(sf_object, feature_name) {
  # Drop the geometry column
  df <- sf::st_drop_geometry(sf_object)

  # Check if the feature exists in the data frame
  feature_exists <- feature_name %in% base::names(df)

  return(feature_exists)
}


#' @name getPanel
#'
#' @title Extract all points with a specified panel id from a sample.
#'
#' @description This is the main function for selecting sites using the BAS master
#' sample. It assumes that you have already defined the master sample using the
#' BoundingBox() function or will be selecting a marine master sample site in BC.
#'
#' @author Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) containing a sample that
#' contains a feature column named panel_id.
#' @param panelid The overlapped panel in the shapefile shp the user wants
#' sample points from.
#'
#' @return The sample for the specified panel.

#' @examples
#' # Halton frame overlapping panel design showing use of getPanel.
#'
#' # Use the North Carolina shapefile supplied in the sf R package.
#' shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
#' shp_gates <- shp_file[shp_file$NAME == "Gates",]
#'
#' # Vertically aligned master sample bounding box.
#' bb <- spbal::BoundingBox(shapefile = shp_gates)
#'
#' # Three panels, of 20 samples each.
#' panels <- c(20, 20, 20)
#'
#' # second panel overlaps first panel by 5, and third panel
#' # overlaps second panel by 5.
#' panel_overlap <- c(0, 5, 5)
#'
#' # generate the sample.
#' samp <- spbal::HaltonFrame(J = c(4, 3),
#'                            boundingbox = bb,
#'                            panels = panels,
#'                            panel_overlap = panel_overlap,
#'                            shapefile = shp_gates)
#'
#' # get halton frame data from our sample.
#' samp3 <- samp$hf.pts.shp
#' samp3
#'
#' panelid <- 1
#' olPanel_1 <- spbal::getPanel(samp3, panelid)
#'
#' @export
getPanel <- function(shapefile, panelid){
  # validate our parameters. Ensure panelid is numeric and has a value greater than zero.
  validate_parameters("panelid", c(panelid))

  # Usage: Check if the sf object contains a feature named "panel_id"
  feature_name <- "panel_id"
  exists <- contains_feature(shapefile, feature_name)
  if(!exists){
    base::stop("spbal(getPanel) Simple file object does not contain a feature named panel_id.")
  }

  indx <- base::c()
  for(k in 1:base::length(shapefile$panel_id)){
    tmp <- shapefile$panel_id[k]
    if(base::any(panelid %in% base::unlist(tmp))){
      indx <- base::c(indx, k)
    }
  }
  smp <- shapefile[indx,]

  # assign the spbal attribute to the sample being returned, i.e. the function that created it.
  base::attr(smp, "spbal") <- "getPanel"

  result <- base::list(sample  = smp,
                       indices = indx)
  return(result)
}
