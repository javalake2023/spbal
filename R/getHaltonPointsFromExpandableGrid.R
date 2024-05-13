# getHaltonPointsFromExpandableGrid.R

#' @name getHaltonPointsFromExpandableGrid
#'
#' @title Generate a Halton frame.
#'
#' @description Find the requested number of Halton points from within a study area using the
#' supplied J and seeds parameters. If the number of points are not found on the first attempt,
#' the frame is expanded, and spbal::getHaltonFrame is called again. This process is repeated
#' until the requested number of points are found. The points and the seeds used to generate
#' the sample are returned to the caller.
#'
#' @author Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) of the study area(s).
#' @param N Number of sites to select. If using stratification it is a named vector containing
#' sample sizes of each group.
#' @param J The number of grid cells. A list of 2 values. The default value is c(3, 2).
#' @param bases Co-prime base for the Halton Sequence. The default value is c(2, 3).
#' @param seeds A vector of 2 seeds, u1 and u2.
#' @param crs Coordinate reference system for the shapefile.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#' @param stratify_found_first A flag to indicate whether we have found the first point in the
#' study region or not. Default FALSE.
#'
#' @return A list containing five variables:
#'
#' \itemize{
#' \item \code{i} The index of the first point chosen at random in the study area.
#' \item \code{diff_} Halton points, the intersection of the bounding box and the study area.
#' \item \code{pts.shp} Halton frame, the sample points within the study area.
#' \item \code{bb.new} The bounding box.
#' \item \code{seeds} The u1 and u2 seeds used to generate the sample.
#' }
#'
#' @keywords internal
getHaltonPointsFromExpandableGrid <- function(shapefile,
                                              N,
                                              J = base::c(4, 3),
                                              bases,
                                              seeds,
                                              crs,
                                              verbose = FALSE,
                                              stratify_found_first = FALSE){
  # find first point in study area...
  # Initialise variables.
  i <- 0
  bases <- base::c(2, 3)
  crs <- sf::st_crs(shapefile)
  pts_in_intersection <- 0

  if(!stratify_found_first){
    while(pts_in_intersection < 1){
      # get the frame.
      #message("pt. 1")
      #message(as.numeric(Sys.time())*1000, digits=15)
      result <- getHaltonFrame(shapefile = shapefile,
                               J = J,
                               i = i,
                               bases = bases,
                               seeds = seeds,
                               crs = crs)
      # how many points in the intersection...
      # just get what we need from result.
      hf_ <- result$hf_
      pts.shp <- result$pts.shp
      seeds <- result$seeds

      tmp <- sf::st_cast(pts.shp, "POINT")
      tmp <- sf::st_as_sf(tmp)
      pts <- hf_$halton_frame
      tmp$ID <- base::seq(1, base::dim(pts)[1])

      # find points in the study area.
      # getHaltonFrame will always return with points from the study area.
      #message("pt. 2")
      #message(as.numeric(Sys.time())*1000, digits=15)
      #diff_ <- sf::st_intersection(tmp, shapefile)

      sel_sgbp = st_intersects(x = tmp, y = shapefile)
      sel_logical = lengths(sel_sgbp) > 0
      diff_ = tmp[sel_logical, ]
      #browser()
      #diff_ = shapefile[sel_logical, ]

      #message("pt. 3")
      #message(as.numeric(Sys.time())*1000, digits=15)

      diff_ <- diff_[base::order(diff_$ID),]
      #message("pt. 4")
      #message(as.numeric(Sys.time())*1000, digits=15)

      # expand the Halton frame until we have seeds that give first point in study region.
      # (just in case - generally will find the first point on the first try).
      i <- i + 1
      # find number of points within our study region
      pts_in_intersection <- base::length(diff_$ID)
    } # End pts_in_intersection < 1

    # select a point in the study area at random.
    base::set.seed(seeds[1] + seeds[2])
    k <- base::sample(pts_in_intersection, 1)
    k <- diff_$ID[k]

    if(verbose){
      msg <- "spbal(getHaltonPointsFromExpandableGrid) Random HF point selected: %s."
      msgs <- base::sprintf(msg, k)
      base::message(msgs)
    }

    # index of first point.
    #k <- first.pt$k
    # generate seeds for the remaining points we need.
    seedshift <- base::c(seeds[1] + k - 1, seeds[2] + k - 1)
    if(verbose){
      msg <- "spbal(getHaltonPointsFromExpandableGrid) New seeds for first point u1 = %s, u2 = %s."
      msgs <- base::sprintf(msg, seedshift[1], seedshift[2])
      base::message(msgs)
    }
    seeds <- seedshift
  } # end if !stratify_found_first

  # Now generate our sample with the new seeds...
  # re-initialise variables.
  i <- 0
  pts_in_intersection <- 0

  # run loop until we have sufficient points...
  while (pts_in_intersection < N){
    # get the frame.
    result <- getHaltonFrame(shapefile = shapefile,
                             J = J,
                             i = i,
                             bases = bases,
                             seeds = seeds,
                             crs = crs)
    hf_ <- result$hf_
    #diff_ <- result$sample     # will always be NULL here.
    pts.shp <- result$pts.shp
    bb.new <- result$bb.new
    seeds <- result$seeds

    pts <- hf_$halton_frame
    tmp <- sf::st_cast(pts.shp, "POINT")
    tmp <- sf::st_as_sf(tmp)
    tmp$ID <- base::seq(1, base::dim(pts)[1])
    # find the points in our study area.
    #diff_ <- sf::st_intersection(tmp, shapefile)
    sel_sgbp = st_intersects(x = tmp, y = shapefile)
    sel_logical = lengths(sel_sgbp) > 0
    diff_ = tmp[sel_logical, ]
    diff_ <- diff_[base::order(diff_$ID), ]

    # find number of points within our shapefile.
    pts_in_intersection <- base::length(diff_$ID)
    if(verbose){
      msg <- "spbal(getHaltonPointsFromExpandableGrid) Points in intersection: %s."
      msgs <- base::sprintf(msg, pts_in_intersection)
      base::message(msgs)
    }

    # expand the Halton frame until we have seeds that give first point in study region.
    i <- i + 1

  } # end while

  if(verbose){
    # display some statistics before returning results.
    msg <- "spbal(getHaltonPointsFromExpandableGrid) %s samples found in %s iterations, using J1=%s and J2=%s."
    msgs <- base::sprintf(msg, pts_in_intersection, i, J[1]+i-1, J[2]+i-1)
    base::message(msgs)
  }

  result <- base::list(i       = i,
                       diff_   = diff_,
                       pts.shp = pts.shp,
                       bb.new  = bb.new,
                       seeds   = seeds)
  return(result)
}

