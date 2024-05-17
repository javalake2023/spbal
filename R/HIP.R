# HIP.R

#' @name hipX1split
#'
#' @title First dimension split.
#'
#' @description Split point pairs using the first dimension.
#'
#' @author Phil Davies, Blair Robertson.
#'
#' @param x1pts First dimension component of point pair.
#' @param HaltonIndex Halton indices for all points in x1Hpts.
#' @param BoxIndex Index of current box to process.
#' @param xlevel The current iteration level.
#' @param x1Hpts First dimension component of Halton point pair.
#'
#' @return A variable called HaltonIndex, the updated Halton indices for all points in x1Hpts.
#'
#' @keywords internal
hipX1split <- function(x1pts, HaltonIndex, BoxIndex, xlevel, x1Hpts) {

  # Determine points in current box
  inBox <- base::which(HaltonIndex == BoxIndex)
  x1pts <- x1pts[inBox]

  # Randomly remove one point (if needed)
  F <- inBox
  nF <- base::length(F)

  # odd number of points?
  if ((nF %% 2) == 1) {
    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x1pts <- x1pts[-r]
  }

  # Partition and update Halton indices
  x1sort <- base::order(x1pts)
  F <- F[x1sort]
  Hpt1 <- x1Hpts[(BoxIndex + 1)]
  Hpt2 <- x1Hpts[(BoxIndex + 1 + xlevel)]
  I <- base::order(base::c(Hpt1, Hpt2))

  # set halton index - top half.
  HaltonIndex[F[1:(base::floor(nF/2))]] <- (BoxIndex + (I[1] - 1) * xlevel)
  # set halton index - bottom half.
  HaltonIndex[F[((base::floor(nF/2))+1):nF]] <- (BoxIndex + (I[2] - 1) * xlevel)
  # return the halton index.
  return(HaltonIndex)
}


#' @name hipX2split
#'
#' @title Second dimension split.
#'
#' @description Split point pairs using the second dimension.
#'
#' @author Phil Davies, Blair Robertson.
#'
#' @param x2pts Second dimension component of point pair.
#' @param HaltonIndex Halton indices for all points in x2Hpts.
#' @param BoxIndex Index of current box to process.
#' @param xlevel The current iteration level.
#' @param x2Hpts Second dimension component of Halton point pair.
#'
#' @return A variable called HaltonIndex, the updated Halton indices for all points in x2Hpts.
#'
#' @keywords internal

hipX2split <- function(x2pts, HaltonIndex, BoxIndex, xlevel, x2Hpts) {

  # Determine points in current box
  inBox <- base::which(HaltonIndex == BoxIndex)
  x2pts <- x2pts[inBox]

  # Randomly remove one or two points (if needed)
  F <- inBox
  # get number of points.
  nF <- base::length(F)

  # not divisible by 3.
  if ((nF %% 3) == 1) {
    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]
  } else if ((nF %% 3) == 2) {
    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]

    r <- base::sample(nF, 1)
    HaltonIndex[F[r]] <- NA
    F <- F[-r]
    nF <- (nF - 1)
    x2pts <- x2pts[-r]
  }

  # Partition and update Halton indices
  x2sort <- base::order(x2pts)
  F <- F[x2sort]
  Hpt1 <- x2Hpts[(BoxIndex + 1)]
  Hpt2 <- x2Hpts[(BoxIndex + 1 + xlevel)]
  Hpt3 <- x2Hpts[(BoxIndex + 1 + (2 * xlevel))]
  I <- base::order(base::c(Hpt1, Hpt2, Hpt3))

  # update halton indices - top third.
  HaltonIndex[F[1:(base::floor(nF/3))]] <- (BoxIndex + (I[1] - 1) * xlevel)
  # update halton indices - middle third.
  HaltonIndex[F[((base::floor(nF/3))+1):(2*(base::floor(nF/3)))]] <- (BoxIndex + (I[2] - 1) * xlevel)
  # update halton indices - bottom third.
  HaltonIndex[F[((2*(base::floor(nF/3)))+1):nF]] <- (BoxIndex + (I[3] - 1) * xlevel)
  # return the halton index.
  return(HaltonIndex)
}


#' @name hipPartition
#'
#' @title Partition the population.
#'
#' @description Partition the resource into boxes with the same nested structure as Halton boxes.
#' The **spbal** parameter **iterations** defines the number of boxes used in the HIP partition and
#' should be larger than the sample size but less than the population size.
#'
#' @author Phil Davies, Blair Robertson.
#'
#' @param pts The population of points.
#' @param its The number of partitioning iterations.
#'
#' @return A list containing the following variables:
#'
#' \itemize{
#' \item \code{ptsIndex} The population index.
#' \item \code{HaltonIndex} Updated Halton indices for all points in pts.
#' }
#'
#' @keywords internal

hipPartition <- function(pts, its) {
  # Initialize
  N <- base::nrow(pts)
  pts <- base::cbind(pts, 1:N)
  # Generate N Halton points.
  res <- cppRSHalton_br(n = N)
  Hpts <- res$pts
  HaltonIndex <- base::rep(0, N)

  # Partitioning parameters
  xlevel <- base::c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184, 10368, 20736)
  #xlevel <- c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184)
  # support upto 13 levels.
  diml <- base::c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1)
  #diml <- c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1)

  # Partitioning loop
  for (j in (1:its)) {
    # for every level...
    if (diml[j] == 1) {
      level <- (xlevel[j] - 1)
      for (i in 0:level) {
        # go calculate halton indices.
        HaltonIndex <- hipX1split(pts[,1], HaltonIndex, i, xlevel[j], Hpts[,1])
      }
    } else {
      level <- (xlevel[j] - 1)
      for (i in 0:level) {
        # go calculate halton indices.
        HaltonIndex <- hipX2split(pts[,2], HaltonIndex, i, xlevel[j], Hpts[,2])
      }
    }

    # Remove discarded points
    TF <- base::which(!base::is.na(HaltonIndex))
    HaltonIndex <- HaltonIndex[TF]
    pts <- pts[TF,]
  }
  ptsIndex <- pts[,3]
  return(base::list(ptsIndex = ptsIndex,
                    HaltonIndex = HaltonIndex))
}


#' @name hipIndexRandomPermutation
#'
#' @title Permute Halton indices.
#'
#' @description A description.
#'
#' @author Phil Davies, Blair Robertson.
#'
#' @param its The number of partitioning iterations.
#'
#' @return A list containing the following variables:
#'
#' \itemize{
#' \item \code{permHaltonIndex} The permuted halton indices for all points.
#' \item \code{B} The number of Halton boxes.
#' }
#'
#' @keywords internal

hipIndexRandomPermutation <- function(its) {

  # Partitioning parameters - support upto 13 levels.
  diml <- base::c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1)
  #diml <- c(1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1)
  diml <- diml[1:its]

  # Halton Indices
  b1 <- 2
  b2 <- 3
  J1 <- base::sum(diml == 1)
  J2 <- base::sum(diml == 2)
  B <- (b1^J1) * (b2^J2)
  HIM <- base::matrix(0, nrow = (b2^J2), ncol = (b1^J1))
  #H <- HaltonPts(B)
  #H <- spbal::cppRSHalton(B)
  #H <- H[,2:3]
  res <- cppRSHalton_br(n = B)
  H <- res$pts

  Hindex <- base::floor(base::rbind((b1^(J1) + 1e-12) * H[,1], (b2^(J2) + 1e-12) * H[,2])) + 1
  Hindex <- base::t(Hindex)

  # Halton Matrix
  for (i in 0:(B - 1)) {
    HIM[((b2^J2) + 1) - Hindex[(i + 1), 2], Hindex[(i + 1), 1]] <- i
  }

  # Permutated Halton Matrix
  step2 <- base::c(2, 4, 8, 16, 32, 64, 128, 256)
  step3 <- base::c(3, 9, 27, 81, 243)
  order2 <- base::c((base::sample(2) - 1), base::rep(NA, ((b1^J1) - b1)))
  order3 <- base::c((base::sample(3) - 1), base::rep(NA, ((b2^J2) - b2)))

  for (i in 1:(J1 - 1)) {
    if (J1 < 2) {
      next
    }
    v <- base::vector()
    L <- base::sum(!base::is.na(order2))

    for (j in (1:L)) {
      k <- order2[j]
      P <- (base::sample(2) - 1)
      s <- step2[i]
      v <- c(v, (k + s * P[1]), (k + s * P[2]))
    }
    order2[1:base::length(v)] <- v
  }
  for (i in 1:(J2 - 1)) {
    if (J2 < 2) {
      next
    }
    v <- base::vector()
    L <- base::sum(!base::is.na(order3))
    for (j in (1:L)) {
      k <- order3[j]
      P <- (base::sample(3) - 1)
      s <- step3[i]
      v <- c(v, (k + s * P[1]),  (k + s * P[2]),  (k + s * P[3]))
    }
    order3[1:base::length(v)] <- v
  }

  # Transform the matrix
  b2vals <- (HIM[1,] %% (b1^J1))
  b3vals <- (HIM[,1] %% (b2^J2))

  I2 <- vector()
  for (i in 1:(b1^J1)) {
    F <- base::which(b2vals == order2[i])
    I2 <- c(I2, F)
  }
  I3 <- vector()
  for (i in 1:(b2^J2)) {
    F <- base::which(b3vals == order3[i])
    I3 <- c(I3, F)
  }

  permHIM <- HIM[I3, I2]
  permHaltonIndex <- base::rep(0, B)

  for (i in 0:(B - 1)) {
    permHaltonIndex[i + 1] <- permHIM[base::which(HIM == i)]
  }

  return(base::list(permHaltonIndex = permHaltonIndex,
                    B               = B))
}


#' @name is_sf_points
#'
#' @title Check if an object is an sf points object.
#'
#' @description Tests if the object passed to the function is a sf points object or not.
#' An internal only function.
#'
#' @author Phil Davies, Blair Robertson.
#'
#' @details Detect if an object is a sf points object or not.
#'
#' @param x A probable sf points object.
#'
#' @return Either TRUE or FALSE.
#'
#' @keywords internal

is_sf_points <- function(x) {
  base::inherits(x, "sf") && base::inherits(x$geometry, "sfc_POINT")
}


#' @name HIP
#'
#' @title Halton Iterative Partitioning (HIP).
#'
#' @description HIP draws spatially balanced samples and over-samples from point
#' resources by partitioning the resource into boxes with the same nested structure
#' as Halton boxes. The **spbal** parameter **iterations** defines the number of
#' boxes used in the HIP partition and should be larger than the sample size but
#' less than the population size. The **iterations parameter** also defines the
#' number of units available in the HIP over-sample, where the over-sample contains
#' one unit from each box in the HIP partition.
#'
#' @author Phil Davies, Blair Robertson.
#'
#' @details Halton iterative partitioning (HIP) extends Basic acceptance
#' sampling (BAS) to point resources. It partitions the resource into $B >= n$
#' boxes that have the same nested structure as in BAS, but different sizes.
#' These boxes are then uniquely numbered using a random-start Halton sequence
#' of length $B$. The HIP sample is obtained by randomly drawing one point from
#' each of the boxes numbered $1, 2, . . . , n$.
#'
#' @param population A population of point pairs.
#' @param n The number of points to draw from the population. Default 20.
#' @param iterations The levels of partitioning required. Default 7.
#' @param minRadius If specified, the minimum distance, in meters, allowed between sample
#' points. This is applied to the $overSample.
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
#' @param verbose Boolean if you want to see any output printed to screen. Helpful if taking a
#' long time. Default is FALSE i.e. no informational messages are displayed.
#'
#' @return Return a list containing the following five variables:
#'
#' \itemize{
#' \item \code{Population} Original population point pairs as an sf object.
#' \item \code{HaltonIndex} The Halton index for the point. Points will be spread equally across
#' all Halton indices.
#' \item \code{sample} The sample from the population of size n.
#' \item \code{overSample} The overSample contains one point from each Halton box. All contiguous
#' sub-samples from oversample are spatially balanced, and the first n points are identical to sample.
#' \item \code{minRadius} This result variable contains the sample created using the minRadius
#' parameter. If the minRadius parameter is not specified then the minRadius variable will contain NULL.
#' }
#'
#' @examples
#' # generating 20 points from a population of 5,000 (random) points with 7
#' # levels of partitioning (4 in the first dimension and 3 in the second) to
#' # give (2^4) * (3^3) = 32 * 27, resulting in 864 boxes ----------------------
#'
#' # set random seed
#' set.seed(511)
#'
#' # define HIP parameters.
#' pop <- matrix(runif(5000*2), nrow = 5000, ncol = 2)
#' n <- 20
#' its <- 7
#'
#' # Convert the population matrix to an sf point object.
#' sf_points <- sf::st_as_sf(data.frame(pop), coords = c("X1", "X2"))
#' dim(sf::st_coordinates(sf_points))
#'
#' # generate HIP sample.
#' result <- spbal::HIP(population = sf_points,
#'                      n = n,
#'                      iterations =  its)
#'
#' # HaltonIndex
#' HaltonIndex <- result$HaltonIndex
#' table(HaltonIndex)
#'
#' # Population Sample
#' HIPsample <- result$sample
#' HIPsample
#'

#' @export
HIP <- function(population = NULL,
                n = 20,
                iterations = 7,
                minRadius = NULL,
                panels = NULL,
                panel_overlap = NULL,
                verbose = FALSE) {

  if(base::is.null(population)){
    base::stop("spbal(HIP) The population parameter must be used. Please specify a population.")
  }
  # must be numeric, greater than 0.
  validate_parameters("hipN", c(n))
  # must be numeric, greater than 1 and less than or equal to 13.
  validate_parameters("hipIterations", c(iterations))
  if(!base::is.null(minRadius)){
    # must be numeric, greater >= 1.
    validate_parameters("minRadius", c(minRadius))
  }

  # population can be either a matrix or sf point object
  sf_points_object <- FALSE
  if (is_sf_points(population)) {
    # It's an sf points object!
    sf_points_object <- TRUE
    # save original population
    sf_population <- population
    # just get coordinates from the population
    population <- sf::st_coordinates(population)
  } else if (base::is.matrix(population)) {
    # It's not an sf points object.
    # must be numeric.
    validate_parameters("hipPopulation", c(population))

    # only 2-d matrix is supported
    if(base::dim(population)[2] != 2){
      stop(base::c("spbal(HIP) Parameter population must have two dimensions."))
    }
  } else {
    # unsupported object type.
    stop(base::c("spbal(HIP) Unsupported type, the population must be defined as a sf point object or matrix."))
  }

  # Partitioning parameters
  iteration_level <- base::c(1, 2, 6, 12, 24, 72, 144, 432, 864, 1728, 5184, 10368, 20736)
  if(((base::length(population)/2)/iteration_level[iterations+1]) < 1.0){
    msg <- "spbal(HIP) Pop. size %s not compatible with number of iteration levels %s. Try pop. size of %s+ or a smaller iteration level."
    msgs <- base::sprintf(msg, base::length(population)/2, iterations, iteration_level[iterations+1])
    stop(msgs)
  }

  # validate panel design if we are using one.
  res <- ValidatePanelDesign(panels, panel_overlap, n)
  panel_design  <- res$panel_design
  number_panels <- res$number_panels
  panel_overlap <- res$panel_overlap
  n             <- res$n

  # Partition the population
  partitionResult <- hipPartition(population, iterations)
  popIndex <- partitionResult$ptsIndex
  HaltonIndex <- partitionResult$HaltonIndex

  # Random permutation of Halton indices
  permResult <- hipIndexRandomPermutation(iterations)
  permHaltonIndex <- permResult$permHaltonIndex
  B <- permResult$B

  # initialise the Order variable.
  Order <- base::rep(0, base::length(HaltonIndex))

  for (i in 0:(B - 1)) {
    Order[base::which(HaltonIndex == i)] <- permHaltonIndex[i + 1]
  }

  # Assign unique indices
  ptsInBox <- (base::length(HaltonIndex) / B)
  if (ptsInBox > 1) {
    for (i in 0:(B - 1)) {
      F <- base::which(Order == i)
      Order[F] <- (i + B * (base::sample(ptsInBox) - 1))
    }
  }

  # Sample Indices
  sampleI <- base::rep(1, B) # was n
  hI <- base::rep(0, B) # was n
  for (i in (1:B)) { # was n
    sampleI[i] <- popIndex[base::which(Order == (i - 1))]
    hI[i] <- (i - 1)
  }

  # add extra column to population sf point object
  sf_points <- NULL
  if(sf_points_object){
    sf_points <- population
  }
  # create PopulationSample from sampleI and return as a sf point object.
  PopulationSample <- sf_population[sampleI,]

  # need to add the following to the sf_population sf points object:
  # the Order and HaltonIndex.
  PopulationSample$HaltonIndex <- hI
  PopulationSample$spbalSeqID <- base::seq(1, B) # was n

  # sort PopulationSample on HaltonIndex
  PopulationSample <- PopulationSample[base::order(PopulationSample$HaltonIndex),]

  # go assign panelid's if required. (returns $sample)
  result <- PanelDesignAssignPanelids(PopulationSample, panels, panel_overlap, panel_design, number_panels)

  # specific column ordering in result$sample.
  fixed_order <- base::c("HaltonIndex", "spbalSeqID")
  remaining_cols <- base::names(result$sample)[-base::match(fixed_order, base::names(result$sample))]
  new_col_order <- base::c(fixed_order, remaining_cols)
  sample <- result$sample[, new_col_order]

  PopulationSample <- sample[1:n,]
  overSample <- sample

  # go filter overSample if user has specified a minimum radius.
  S <- filterOnDistance(overSample, minRadius)

  # assign the spbal attribute to the sample being returned, i.e. the function that created it.
  base::attr(PopulationSample, "spbal") <- "HIP"
  base::attr(overSample, "spbal") <- "HIP"

  return(list(Population       = sf_points,
              HaltonIndex      = HaltonIndex,
              sample           = PopulationSample,
              overSample       = overSample,
              minRadius        = S))
}
