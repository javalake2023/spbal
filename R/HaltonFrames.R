# HaltonFrames.R

#' @import units
#' @import sf
#' @import Rcpp
#' @useDynLib spbal, .registration = TRUE

#' @name HaltonFrame
#'
#' @title Create a Halton Frame.
#'
#' @description Halton frames discretize an areal resource into a spatially ordered grid,
#' where samples of consecutive frame points are spatially balanced. To generate Halton Frames,
#' spbal requires a study region \code{shapefile} and the region’s \code{bounding box}.
#'
#' @author Phil Davies.
#'
#' @param N The number of points in the frame to generate.
#' @param J The number of grid cells. A list of 2 values. The default value is c(3, 2).
#' @param bases Co-prime base for the Halton Sequence. The default value is c(2, 3).
#' @param shapefile A sf object. If the shapefile parameter is NULL then function
#' HaltonFrameBase is called directly.
#' @param boundingbox Bounding box around the study area. If a bounding box is not supplied
#' then spbal will generate a bounding box for the shapefile.
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
#' @param seeds A vector of 2 seeds, u1 and u2. If not specified, the default is NULL.
#' @param stratum Name of column in shapefile that makes up the strata.
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return Returns a list containing five variables:
#'
#' \itemize{
#' \item \code{J} The number of grid cells. A list of 2 values that were used to generate this
#' Halton grid and frame.
#' \item \code{hg.pts.shp} Halton grid over the bounding box and study area.
#' \item \code{hf.pts.shp} Halton frame, the sample points within the study area.
#' \item \code{bb} The bounding box.
#' \item \code{seeds} The u1 and u2 seeds used to generate the sample.
#' }
#'
#' The sample points in \code{hf.pts.shp} are returned in the form of a simple feature
#' collection of POINT objects. As well as having the features from the original \code{shapefile},
#' the following new attributes have been added:
#'
#' \itemize{
#'   \item \code{spbalSeqID}: A unique identifier for every sample point.
#'   \item \code{ID}: A unique identifier, the Halton frame point order.
#' }
#'
#' @examples
#' # we discretize the Gates study region into a coarse grid using
#' # B = 2^{J_1} * 3^{J_2}= (2^3) * (3^2) (9 by 8 grid) ------------------------
#'
#' # Use the North Carolina shapefile supplied in the sf R package.
#' shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
#' shp_gates <- shp_file[shp_file$NAME == "Gates",]
#'
#' # Vertically aligned master sample bounding box.
#' bb <- spbal::BoundingBox(shapefile = shp_gates)
#'
#' set.seed(511)
#' result6 <- spbal::HaltonFrame(shapefile = shp_gates,
#'                               J = c(3, 2),
#'                               boundingbox = bb)
#' # get the frame points.
#' Frame <- result6$hf.pts.shp
#' Frame
#' # get the grid points.
#' Grid <- result6$hg.pts.shp
#' Grid
#'
#' @export
HaltonFrame <- function(N = 1,
                        J = base::c(3, 2),
                        bases = base::c(2, 3),
                        boundingbox = NULL,
                        shapefile = NULL,
                        panels = NULL,
                        panel_overlap = NULL,
                        seeds = NULL,
                        stratum = NULL,
                        verbose = FALSE){

  # initialize variables.
  # defaults, before we see what the user wants.
  wantHaltonGrid <- FALSE
  wantHaltonFrame <- FALSE
  hf_stratification <- FALSE

  if(base::is.null(N)){
    wantHaltonGrid <- TRUE
    wantHaltonFrame <- FALSE
    if(verbose){
      base::message("spbal(HaltonFrame) Request for a Halton Grid.")
    }
  }

  # validate the shapefile parameter.
  if(base::is.null(shapefile)){
    base::stop("spbal(HaltonFrame) The shapefile parameter must be used. Please specify a shapefile.")
  }

  # validate other parameters.
  validate_parameters("J", J)
  validate_parameters("bases", bases)
  if (!base::is.null(N)){
    validate_parameters("N", base::c(N))
  }

  # A bounding box must be specified.
  # If user has not specified a bb then we will generate one.
  if(base::is.null(boundingbox)){
    boundingbox <- spbal::BoundingBox(shapefile = shapefile)
  }

  if (base::is.null(seeds)){
    seeds <- generateUVector()
  } else {
    validate_parameters("seeds", seeds)
  }

  # validate panel design if we are using one.
  res <- ValidatePanelDesign(panels, panel_overlap, N)
  panel_design  <- res$panel_design
  number_panels <- res$number_panels
  panel_overlap <- res$panel_overlap
  N             <- res$n

  # if both not NULL then we want stratification.
  if(!base::is.null(base::names(N)) & !base::is.null(stratum)){
    hf_stratification <- TRUE
    wantHaltonFrame <- TRUE
    strata.levels <- base::names(N)
    if(verbose){
      msg <- "spbal(HaltonFrame) Stratification request for the following strata: %s.\n"
      msgs <- base::sprintf(msg, strata.levels)
      base::message(msgs)
    }
  } else {
    if(N >= 1){
      wantHaltonGrid <- TRUE
      wantHaltonFrame <- TRUE
      if(verbose){
        # state how many samples user is looking for.
        msg <- "spbal(HaltonFrame) Request for %s samples from a Halton Frame."
        msgs <- base::sprintf(msg, N)
        base::message(msgs)
      }
    }
  }

  # ensure the shapefile (if specified) has an associated CRS.
  crs <- NULL
  if (!base::is.null(shapefile)){
    if (base::is.null(sf::st_crs(shapefile))) {
      stop("spbal(HaltonFrame) Shapefile does not have an associated CRS.")
    } else {
      if(verbose){
        msg <- "spbal(HaltonFrame) Shapefile has an associated CRS."
        msgs <- base::sprintf(msg)
        base::message(msgs)
      }
      crs <- sf::st_crs(shapefile)
    }
  } else {
    # shapefile is null, ie. not specified, so just call HaltonFrameBase
    hf_ <- HaltonFrameBase(J = base::c(J[1], J[2]), bases = bases, seeds = seeds)
    return(hf_)
  }

  # number of points currently in the area of interest.
  pts_in_intersection <- 0
  # initialise
  i <- 0

  if(wantHaltonGrid & !wantHaltonFrame){
    # go get halton grid.
    result <- getHaltonFrame(shapefile, J, i, bases, seeds, crs)
    hf_ <- result$hf_
    diff_ <- result$sample
    pts.shp <- result$pts.shp
    bb.new <- result$bb.new
    seeds <- result$seeds

  } else {

    # need to check hf_stratification before we run the while loop.
    # run the while loop in new function so we can loop for each strata level for the desired N.
    # save the returned points using rbind.
    if(hf_stratification){
      # perform stratification
      smp <- NULL
      stratification_found_first_point <- FALSE

      for(h in 1:base::length(N)){
        if(verbose){
          msg <- "spbal(HaltonFrame) Stratum: %s."
          msgs <- base::sprintf(msg, strata.levels[h])
          base::message(msgs)
        }

        h.indx <- base::which(shapefile[, stratum, drop = TRUE] == strata.levels[h])
        shp.stratum <- shapefile[h.indx,]

        # find the first point in the study region (picked at random)
        #first.pt <- findFirstStudyRegionHFPoint(shapefile = shp.stratum, bb = boundingbox, seeds = seeds, verbose = verbose)
        # index of first point.
        #k <- first.pt$k
        # generate seeds for the remaining points we need.
        #seedshift <- base::c(first.pt$seeds[1] + k - 1, first.pt$seeds[2] + k - 1)

        result <- getHaltonPointsFromExpandableGrid(shapefile = shp.stratum,
                                                    N = N[h],
                                                    J = J,
                                                    bases = bases,
                                                    seeds = seeds,
                                                    crs = crs,
                                                    verbose = verbose,
                                                    stratify_found_first = stratification_found_first_point)

        stratification_found_first_point <- TRUE
        #seedshift <- result$seed
        seeds <- result$seed
        diff_pts <- result$diff_
        df_sorted <- diff_pts[base::order(diff_pts$ID),]

        #ret_sample <- base::rbind(first.pt$first.pt, df_sorted)
        #sorted_samp <- ret_sample
        sorted_samp <- df_sorted
        sorted_samp$spbalSeqID <- base::seq(1, base::length(sorted_samp$ID))
        # add stratum as a new column.
        sorted_samp[[stratum]] <- strata.levels[h]

        diff_ <- sorted_samp[1:N[h],]
        # return original seeds.
        #seeds <- first.pt$seeds
        #
        #df_sorted$spbalSeqID <- base::seq(1, base::length(df_sorted$ID))
        # return N[h] from the intersection for current stratum.
        #diff_ <- df_sorted[1:N[h],]
        smp <- base::rbind(smp, diff_)

      } # end for h

      # load variables for return to caller.
      i <- result$i
      diff_ <- smp
      pts.shp <- result$pts.shp
      bb.new <- result$bb.new
      seeds <- result$seeds
      #seeds <- first.pt$seeds

    } else {
      # stratification not required - use entire study area.

      # find the first point in the study region (picked at random)
      #first.pt <- findFirstStudyRegionHFPoint(shapefile = shapefile, bb = boundingbox, seeds = seeds, verbose = verbose)
      # index of first point.
      #k <- first.pt$k
      # generate seeds for the remaining points we need.
      #seedshift <- base::c(first.pt$seeds[1] + k - 1, first.pt$seeds[2] + k - 1)
      # go get the Halton points.
      result <- getHaltonPointsFromExpandableGrid(shapefile = shapefile,
                                                  N = N,
                                                  J = J,
                                                  bases = bases,
                                                  seeds = seeds,
                                                  crs = crs,
                                                  verbose = verbose,
                                                  stratify_found_first = FALSE)
      i         <- result$i
      diff_     <- result$diff_
      pts.shp   <- result$pts.shp
      bb.new    <- result$bb.new
      seedshift <- seeds         # was seeds until first.pt code was added. first.pt$seeds
    }

  } # end if (wantHaltonGrid & !wantHaltonFrame

  # are we performing a panel_design? yes then go assign panelid's.
  if(panel_design){
    if(verbose){
      base::message("spbal(HaltonFrame) User requested a Panel Design.")
    }
    diff_pts <- sf::st_cast(diff_, "POINT")
    df_sorted <- diff_pts[base::order(diff_pts$ID), ]
    df_sorted$spbalSeqID <- base::seq(1, base::length(df_sorted$ID))
    diff_pts_sf <- sf::st_as_sf(df_sorted)
    res <- PanelDesignAssignPanelids(diff_pts_sf, panels, panel_overlap, panel_design, number_panels)
    diff_ <- res$sample
  }

  # if we are not performing a panel_design or stratification
  if (!panel_design & wantHaltonFrame & !hf_stratification){
    if(verbose){
      base::message("spbal(HaltonFrame) Return N sample points.")
    }
    # turn our sample into points.
    diff_pts <- sf::st_cast(diff_, "POINT")
    df_sorted <- diff_pts[base::order(diff_pts$ID), ]
    # return everything from the intersection.
    diff_ <- df_sorted #[1:n,]

    # handle the Halton grid.
    hg_pts <- sf::st_cast(pts.shp, "POINT")
    hg_pts <- sf::st_as_sf(hg_pts)
    hg_pts$ID <- seq(1, length(hg_pts$x))
    pts.shp <- hg_pts

    # handle first.pt
    #f.pt <- sf::st_cast(first.pt$first.pt, "POINT")
    #f.pt <- sf::st_as_sf(f.pt)
    #zzz <- sf::st_as_sf(base::data.frame(SiteID = f.pt$ID, f.pt$x))
    #ret_sample <- base::rbind(first.pt$first.pt, diff_)
    #sorted_samp <- ret_sample[base::order(ret_sample$SiteID), ]
    #sorted_samp <- ret_sample
    sorted_samp <- diff_
    sorted_samp$spbalSeqID <- base::seq(1, base::length(sorted_samp$ID))
    if(N == 1){
      diff_ <- sorted_samp #[1:N,]
    } else {
      diff_ <- sorted_samp[1:N,]
    }
    # return original seeds.
    seeds <- seeds # first.pt$seeds
  }

  # specific column ordering in diff_.
  fixed_order <- base::c("ID", "spbalSeqID")
  # re-order columns for the frame.
  remaining_cols <- base::names(diff_)[-base::match(fixed_order, base::names(diff_))]
  new_col_order <- base::c(fixed_order, remaining_cols)
  diff_ <- diff_[, new_col_order]
  # no need to re-order on the grid as we display all the points, not just the study area.

  # assign the spbal attribute to the sample being returned, i.e. the function that created it.
  base::attr(diff_, "spbal") <- "HaltonFrame"
  base::attr(pts.shp, "spbal") <- "HaltonFrame"

  # Need to return cpprshs$pts, cpprshs$xklist, z and hf
  result <- base::list(J          = c(J[1]+i-1, J[2]+i-1),
                       hf.pts.shp = diff_,   # Halton Frame
                       hg.pts.shp = pts.shp, # Halton Grid
                       bb         = bb.new,
                       seeds      = seeds)
  return(result)
}


#' @name HaltonFrameBase
#'
#' @title Generate a Halton Frame.
#'
#' @description A description of this useful function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param n The number of points in the frame to generate.
#' @param J The number of grid cells. A list of 2 values. The default value is c(3, 2), we could also use c(5, 3).
#' @param bases Co-prime base for the Halton Sequence. The default value is c(2, 3).
#' @param seeds The u1 and u2 seeds to use.
#'
#' @return A list containing the following four variables:
#' halton_seq -
#' halton_seq_div -
#' Z -
#' halton_frame -
#'
#' @keywords internal
HaltonFrameBase <- function(n = (bases[1]^J[1]) * (bases[2]^J[2]),
                            J = base::c(3, 2),
                            bases = base::c(2, 3),
                            seeds = NULL){
  # validate our parameters.
  #validate_parameters("J", J)
  #validate_parameters("bases", bases)
  #validate_parameters("n", base::c(n))
  #if (!base::is.null(seeds)){
  #  validate_parameters("seeds", seeds)
  #}

  #
  j1 <- J[1]
  j2 <- J[2]
  # calculate B
  B <- (bases[1]^j1) * (bases[2]^j2)
  # check how many points the caller wants.
  if(n > B) {
    B <- (base::floor(n / B) + 1) * B
  }
  # double the number of points
  B2 <- 2 * B
  # compute B2 halton points
  #cpprshs <- spbal::cppRSHalton_br(n = B2, bases = bases)
  if(base::is.null(seeds)){
    #cpprshs <- cppBASpts(n = B2, bases = bases)
    cpprshs <- cppRSHalton_br(n = B2, bases = bases)
  } else {
    #cpprshs <- cppBASpts(n = B2, bases = bases, seeds = seeds)
    cpprshs <- cppRSHalton_br(n = B2, bases = bases, seeds = seeds)
  }
  # pull points closer to center.
  z <- ((cpprshs$pts[1:B,] + cpprshs$pts[(B+1):B2,])/2)
  # x-dimension
  x_dim <- (base::floor(bases[1]^j1 * z[,1])/bases[1]^j1) + 0.5*(1/(bases[1]^j1))
  # y-dimension
  y_dim <- (base::floor(bases[2]^j2 * z[,2])/bases[2]^j2) + 0.5*(1/bases[2]^j2)
  # Halton Frame
  hf <- base::cbind(x_dim, y_dim)
  # Create an index number for each data point.
  halton_indx <- base::seq(1, B2/2)

  # Need to return cpprshs$pts, cpprshs$xklist, z, hf and seeds.
  result <- base::list(halton_seq     = cpprshs$pts,
                       halton_seq_div = cpprshs$xklist,
                       Z              = z,
                       halton_frame   = hf,
                       halton_indx    = halton_indx,
                       seeds          = cpprshs$seeds)
  return(result)
}


#' @name getHaltonFrame
#'
#' @title Obtain a Halton Frame over a shapefile.
#'
#' @description An internal only function.
#'
#' @details This function was written by Phil Davies.
#'
#' @param shapefile A MULTIPOINT or POINT object that we want to generate a halton frame for.
#' @param J The number of grid cells. A list of 2 values.
#' @param bases Co-prime base for the Halton Sequence.
#' @param i An integer to add to the J parameter elements to expand the Halton Frame in both
#' directions if the required number of sample points cannot be found in the region of interest
#' in the current Halton frame.
#' @param seeds A list of 2 seeds, u1 and u2.
#' @param crs Coordinate reference system for the shapefile.
#'
#' @return A list containing the following variables: hf_, sample, pts.shp, bb.new, seeds
#'
#' @keywords internal
getHaltonFrame <- function(shapefile, J, i, bases, seeds, crs){
  #
  #seeds <- c(seeds[1]+i, seeds[2]+i)
  hf_ <- HaltonFrameBase(J = base::c(J[1]+i, J[2]+i), bases = bases, seeds = seeds)

  # process points returned.
  pts <- hf_$halton_frame
  pts <- base::cbind(base::seq(1, base::dim(pts)[1]), pts)

  # save returned seeds in case they have changed (would only change if initially NULL).
  seeds <- hf_$seeds

  #
  bb <- sf::st_as_sfc(sf::st_bbox(shapefile))
  cntrd <- sf::st_centroid(bb)
  bb.rot <- (bb - cntrd) * rot(0) + cntrd
  bb.new <- sf::st_as_sfc(sf::st_bbox(bb.rot))

  #
  base::attr(bb.new, "rotation") <- 0
  base::attr(bb.new, "centroid") <- sf::st_coordinates(cntrd)
  pts.shp <- rotate.scale.coords(coords = pts, bb = bb.new)
  # make sure our shapefile has a CRS (needed for plotting later on).
  sf::st_crs(pts.shp) <- crs
  # always return NULL when just generating a Halton Frame.
  diff_ <- NULL

  result <- base::list(hf_     = hf_,
                       sample  = diff_,
                       pts.shp = pts.shp,
                       bb.new  = bb.new,
                       seeds   = seeds)
  return(result)
}

