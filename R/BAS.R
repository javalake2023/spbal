# BAS.R

#' @name BAS
#'
#' @title Balanced Acceptance Sampling (BAS).
#'
#' @description BAS draws spatially balanced samples from areal resources. To draw BAS samples,
#' spbal requires a study region shapefile and the region’s bounding box. An initial sample size
#' is also needed, which can be easily increased or decreased within spbal for master sampling
#' applications
#'
#' @details This function was first written by Paul van Dam-Bates for the
#' package BASMasterSample and later simplified by Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) to select sites for.
#' @param n Number of sites to select. If using stratification it is a named vector containing
#' sample sizes of each group.
#' @param boundingbox Bounding box which defines the Master Sample. A bounding box must be
#' supplied.
#' @param minRadius If specified, the minimum distance, in meters, allowed between sample
#' points. This is applied to the $sample points. Points that meet the minRadius criteria
#' are retuned in the minRadius output variable.
#' @param panels A list of integers that define the size of each panel in a
#' non-overlapping panels design. The length of the list determines the number of
#' panels required. The sum of the integers in the panels parameter will determine
#' the total number of samples selected, n. The default value for panels is NULL,
#' this indicates that a non-overlapping panel design is not wanted.
#' @param panel_overlap A list of integers that define the overlap into the previous
#' panel. Is only used when the panels parameter is not NULL. The default value for
#' panel_overlap is NULL. The length of panel_overlap must be equal to the length
#' of panels. The first value is always forced to zero as the first panel never
#' overlaps any region.
#' @param stratum The name of a column in the data.frame attached to shapefile that defines
#' the strata of interest.
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL and will
#' be defined randomly using function \code{generateUVector}.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing three variables, \code{$sample} containing locations in the BAS sample,
#' in BAS order, \code{$seeds}, the u1 and u2 seeds used to generate the sample and \code{$minRadius}
#' containing points from $sample that meet the minRadius criteria. If the minRadius
#' parameter is NULL then the $minRadius returned will also be NULL.
#'
#' The sample points are returned in the form of a simple feature collection of POINT objects.
#' They have the following attributes:
#' \itemize{
#'   \item \code{SiteID} A unique identifier for every sample point. This
#'   encodes the BAS order.
#'   \item \code{spbalSeqID} A unique identifier for every sample point. This
#'   encodes the BAS sample order.
#'   \item \code{geometry} The XY co-ordinates of the sample point in the CRS of the original
#'   shapefile.
#' }

#' @examples
#' # Equal probability BAS sample ----------------------------------------------
#'
#' # Use the North Carolina shapefile supplied in the sf R package.
#' shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
#' shp_gates <- shp_file[shp_file$NAME == "Gates",]
#'
#' # Vertically aligned master sample bounding box.
#' bb <- spbal::BoundingBox(shapefile = shp_gates)
#'
#' set.seed(511)
#' n_samples <- 20
#' # Equal probability BAS sample.
#' result <- spbal::BAS(shapefile = shp_gates,
#'                      n = n_samples,
#'                      boundingbox = bb)
#' BAS20 <- result$sample
#' # display first three sample points.
#' BAS20[1:3,]
#'
#' # Increase the BAS sample size ----------------------------------------------
#' n_samples <- 50
#' result2 <- spbal::BAS(shapefile = shp_gates,
#'                       n = n_samples,
#'                       boundingbox = bb,
#'                       seeds = result$seed)
#' BAS50 <- result2$sample
#' BAS50[1:3,]
#'
#' # Check, first n_samples points in both samples must be the same.
#' all.equal(BAS20$geometry, BAS50$geometry[1:20])
#'
#' @export
BAS <- function(shapefile = NULL,
                n = 100,
                boundingbox = NULL,
                minRadius = NULL,
                panels = NULL,
                panel_overlap = NULL,
                stratum = NULL,
                seeds = NULL,
                verbose = FALSE){

  # validate shapefile and other BAS parameters.
  # validate the shapefile parameter.
  if(base::is.null(shapefile)){
    base::stop("spbal(BAS) The shapefile parameter must be used. Please specify a shapefile.")
  }
  shp_geometry <- sf::st_geometry_type(shapefile)
  if (!base::all(shp_geometry %in% c("MULTIPOLYGON", "POLYGON"))){
    msg <- "spbal(BAS) Unsupported geometry in shapefile, %s."
    msgs <- base::sprintf(msg, shp_geometry)
    base::stop(msgs)
  }

  # A bounding box must be specified.
  if(base::is.null(boundingbox)){
    msg <- "spbal(BAS) A Bounding Box must be specified. Use spbal::BoundingBox first."
    msgs <- base::sprintf(msg)
    base::stop(msgs)
  }

  if(!base::is.null(minRadius)){
    # must be numeric, greater >= 1.
    validate_parameters("minRadius", c(minRadius))
  }

  # validate panel design if we are using one.
  res <- ValidatePanelDesign(panels, panel_overlap, n)
  panel_design  <- res$panel_design
  number_panels <- res$number_panels
  panel_overlap <- res$panel_overlap
  n             <- res$n

  # stratification wanted?
  if(base::is.null(stratum)){
    # no stratification, just a simple sample wanted.
    result <- getBASSampleDriver(shapefile = shapefile,
                                 n = n,
                                 bb = boundingbox,
                                 seeds = seeds,
                                 verbose = verbose)
    smp <- result$sample
    seeds <- result$seed
  }else{
    if(base::is.null(base::names(n))) base::stop("spbal(BAS) Need design sample size as n = named vector")

    strata.levels <- base::names(n)
    smp <- NULL

    # for every stratum...
    for(k in 1:base::length(n)){
      if(verbose){
        msg <- "spbal(BAS) Stratum: %s."
        msgs <- base::sprintf(msg, strata.levels[k])
        base::message(msgs)
      }
      k.indx <- base::which(shapefile[, stratum, drop = TRUE] == strata.levels[k])
      shp.stratum <- shapefile[k.indx,]
      result <- getBASSampleDriver(shapefile = shp.stratum,
                                   n = n[k],
                                   bb = boundingbox,
                                   seeds = seeds,
                                   verbose = verbose)
      smp.s <- result$sample
      seeds <- result$seed
      smp.s[stratum] <- strata.levels[k]
      smp <- base::rbind(smp, smp.s)
    } # end for k.
  } # end is.null(stratum).

  # reorder columns.
  #col_order <- base::c("SiteID", "spbalSeqID", "geometry")
  #smp <- smp[, col_order]

  # go assign panelid's if required.
  res <- PanelDesignAssignPanelids(smp, panels, panel_overlap, panel_design, number_panels)

  # go filter res$sample if user has specified a minimum radius.
  S <- NULL
  if(!is.null(minRadius)){
    S <- filterOnDistance(res$sample, minRadius)
  }

  # specific column ordering in res$sample.
  fixed_order <- base::c("SiteID", "spbalSeqID")
  remaining_cols <- base::names(res$sample)[-base::match(fixed_order, base::names(res$sample))]
  new_col_order <- base::c(fixed_order, remaining_cols)
  sample <- res$sample[, new_col_order]

  # assign the spbal attribute to the sample being returned, i.e. the function that created it.
  base::attr(sample, "spbal") <- "BAS"

  # return the sample and the u1, u2 seeds used.
  result <- base::list(sample    = sample,
                       seed      = seeds,
                       minRadius = S)
  return(result)
}


#' @name getBASSampleDriver
#'
#' @title Title.
#'
#' @description This function repeatedly calls function spbal::getBASSample to generate the BAS
#' sample. Once the requested number of points within the intersection of the shapefile and the
#' study area have been obtained, the sample and seeds are returned to the caller.
#'
#' @details This function was written by Phil Davies based on origin code by Paul van Dam-Bates
#' from the BASMasterSample package.
#'
#' @param shapefile Shape file as a polygon (sp or sf) to select sites for.
#' @param bb Bounding box which defines the area around the study area. A bounding box must be
#' supplied.
#' @param n Number of sites to select. If using stratification it is a named vector containing
#' sample sizes of each group.
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL and will
#' be defined randomly.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return A list containing two variables, \code{$sample} containing locations in the BAS sample,
#' in BAS order and \code{$seeds}, the u1 and u2 seeds used to generate the sample.
#'
#' @keywords internal
getBASSampleDriver <- function(shapefile, bb, n, seeds, verbose = FALSE){

  # The assumption is that if seeds are not null then user has previously run this to get
  # a sample and associated seeds.
  if(base::is.null(seeds)){
    seeds <- generateUVector()
    if(verbose){
      msg <- "spbal(getBASSampleDriver) Seeds from generateUVector() u1 = %s, u2 = %s."
      msgs <- base::sprintf(msg, seeds[1], seeds[2])
      base::message(msgs)
    }
    # find the first point in the study region (picked at random), specifically we want to
    # find the seeds that give us the first point in the study region.
    first.pt <- findFirstStudyRegionPoint(shapefile = shapefile, bb = bb, seeds = seeds, verbose = verbose)

    # get the index of the first point (actually the SiteID).
    k <- first.pt$k
    # calculate new seeds.
    seedshift <- base::c(first.pt$seeds[1] + k - 1, first.pt$seeds[2] + k - 1)
    if(verbose){
      msg <- "spbal(getBASSampleDriver) New seeds for first point u1 = %s, u2 = %s."
      msgs <- base::sprintf(msg, seedshift[1], seedshift[2])
      base::message(msgs)
    }
  } else {
    seedshift <- seeds
  } # end is.null(seeds)

  # number of samples required.
  draw <- n * 4
  # just the first point so far, need n.
  num_samples <- 1
  # count number of times we call spbal::getBASSample.
  call.getBASSample.cnt <- 0
  # keep generating BAS samples until we find n sample points in the study area.
  while(num_samples < n){
    # double the number of points to find to try and reduce number of loops.
    draw <- draw * 2
    # go get sample.
    pts.sample <- getBASSample(shapefile = shapefile, bb = bb , n = draw, seeds = seedshift)
    # get sample points.
    ret_sample <- pts.sample$sample
    # how many samples do we have?
    n_samples <- base::length(ret_sample$SiteID)

    if(verbose){
      msg <- "spbal(getBASSampleDriver) after getBASSample n_samples = %s. num_samples = %s"
      msgs <- base::sprintf(msg, n_samples, num_samples)
      base::message(msgs)
    }
    call.getBASSample.cnt <- call.getBASSample.cnt + 1
    num_samples <- n_samples
  } # end while num_samples < n

  if(verbose){
    msg <- "spbal(getBASSampleDriver) Needed %s call(s) to obtain %s samples."
    msgs <- base::sprintf(msg, call.getBASSample.cnt, n)
    base::message(msgs)
  }

  seeds <- seedshift
  # add a sample ID for the user.
  ret_sample$spbalSeqID <- base::seq(1, base::length(ret_sample$SiteID))
  # add the "Count" feature for NZ_DOC (acts as a unique ID for backward compatibility purposes).
  ret_sample$Count <- ret_sample$SiteID

  # return sample and seeds to caller.
  result <- base::list(sample = ret_sample[1:n,],
                       seeds  = seeds)
  return(result)
}


#' @name getBASSample
#'
#' @title Generate the BAS sample.
#'
#' @description This function is repeatedly called from function spbal::getBASSampleDriver
#' to generate a BAS sample.
#'
#' @details This function was written by Phil Davies.
#'
#' @param shapefile Shape file as a polygon (sp or sf) to select sites for.
#' @param bb Bounding box which defines the area around the study area. A bounding box must be
#' supplied.
#' @param n Number of sites to select. If using stratification it is a named vector containing
#' sample sizes of each group.
#' @param seeds A vector of 2 seeds, u1 and u2. seeds must have a value when this function is called.
#'
#' @return A list containing two variables, \code{$sample} containing locations in the BAS sample,
#' in BAS order and \code{$seeds}, the u1 and u2 seeds used to generate the sample.
#'
#' @keywords internal
getBASSample <- function(shapefile, bb, n, seeds){

  if(base::is.null(seeds)){
    msg <- "spbal(getBASSample) The seeds parameter must not be NULL."
    msgs <- base::sprintf(msg)
    base::stop(msgs)
  }
  seedshift <- seeds

  # Scale and shift Halton to fit into bounding box
  bb.bounds <- sf::st_bbox(bb)
  scale.bas <- bb.bounds[3:4] - bb.bounds[1:2]
  shift.bas <- bb.bounds[1:2]

  # go generate n Halton points.
  #res <- cppRSHalton_br(n = n, seeds = seedshift, bases = base::c(2, 3))
  res <- cppBASpts(n = n, seeds = seedshift, bases = base::c(2, 3))
  pts <- res$pts
  siteid <- base::seq(from = 1, to = n, by = 1)
  pts <- base::cbind(siteid, pts)

  xy <- base::cbind(pts[,2]*scale.bas[1] + shift.bas[1], pts[,3]*scale.bas[2] + shift.bas[2])

  pts.coord <- sf::st_as_sf(base::data.frame(SiteID = pts[,1], xy), coords = c(2, 3))

  sf::st_crs(pts.coord) <- sf::st_crs(bb)
  # find the intersection. Generates the same as sf::st_intersection(pts.coord, shapefile)
  pts.intersect <- pts.coord[shapefile,]

  # return the point to the driver.
  result <- base::list(sample = pts.intersect,
                       seeds  = seedshift)
  return(result)
}

