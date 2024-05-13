# ValidatePanelDesign.R

#' @name ValidatePanelDesign
#'
#' @title Validate the panels and panel_overlap parameters.
#'
#' @description This function is used to validate the panels and panel_overlap
#' parameters. The panel_design flag is set TRUE when the panels and/or panel_overlap
#' parameters are not NULL. This is an internal only function.
#'
#' @author Phil Davies.
#'
#' @param panels A list of integers that define the size of each panel in a
#' non-overlapping panels design. The length of the list determines the number of
#' panels required. The sum of the integers in the panels parameter will determine
#' the total number of samples selected, n. The default value for panels is NULL,
#' this indicates that a non-overlapping panel design is not wanted.
#' @param panel_overlap A list of integers that define the overlap into the previous
#' panel. It is only used when the panels parameter is not NULL. The default value for
#' panel_overlap is NULL. The length of panel_overlap must be equal to the length
#' of panels. The first value is always forced to zero as the first panel never
#' overlaps any region.
#' @param n The number of samples required. Only used when panels and panel_overlap are NULL.
#'
#' @return A list containing four variables, they are detailed below.
#'
#' \itemize{
#'   \item \code{n} When the \code{panels} parameter is not null, the \code{n} parameter is set
#'   using the sum of all the panel sizes in \code{panels}.
#'   \item \code{panel_design} A boolean, TRUE, indicates that the user wants a panels design.
#'   \item \code{number_panels} The number of panels specified in the panel design.
#'   \item \code{panel_overlap} Updated panel_overlap vector, the first element is always forced
#'   to zero irrespective of what the user specified.
#' }
#'
#' @keywords internal

ValidatePanelDesign <- function(panels, panel_overlap, n){

  # the default is that it's not a panel design. Will be set to true if either of the
  # panels or panel_overlap parameters are not NULL.
  panel_design <- FALSE
  number_panels <- 0

  # verify panels parameter, must be a list of numerics (if not null).
  # if panels not NULL, then we will ignore the n parameter.
  if(!base::is.null(panels)){
    validate_parameters("panels", panels)
    n <- base::sum(panels)
    panel_design <- TRUE
    number_panels <- base::length(panels)
  }

  # verify panels_overlap parameter, must be a list of numerics (if not null).
  if(!base::is.null(panel_overlap)){
    validate_parameters("panel_overlap", panel_overlap)
    if(!panel_design){
      base::stop("spbal(ValidatePanelDesign) panels parameter must be specified when panel_overlap specified.")
    }
    if(base::length(panels) != base::length(panel_overlap)){
      msg <- "spbal(ValidatePanelDesign) length of panels [%s] must match length of panel_overlap [%s]."
      msgs <- base::sprintf(msg, base::length(panels), base::length(panel_overlap))
      base::stop(msgs)
    }
    panel_design <- TRUE
    # force zero for panel 1.
    panel_overlap[1] <- 0
  }
  result <- base::list(panel_design  = panel_design,
                       number_panels = number_panels,
                       panel_overlap = panel_overlap,
                       n             = n)
  return(result)
}


#' @name PanelDesignAssignPanelids
#'
#' @title Assign panel ids to the samples.
#'
#' @description This function assigns panel id's to each sample based on values in the
#' panels and panel_overlap parameters. This is an internal only function.
#'
#' @author Phil Davies.
#'
#' @param smp The shapefile for the region under study.
#' @param panels A list of integers that defines the size of each panel in a
#' non-overlapping panels design. The length of the list determines the number of
#' panels required. The sum of the integers in the panels parameter will determine
#' the total number of samples selected, n. The default value for panels is NULL,
#' this indicates that a non-overlapping panel design is not wanted.
#' @param panel_overlap A list of integers that define the overlap into the previous
#' panel. Is only used when the panels parameter is not NULL. The default value for
#' panel_overlap is NULL. The length of panel_overlap must be equal to the length
#' of panels. The first value is always forced to zero as the first panel never
#' overlaps any region.
#' @param panel_design A flag, when TRUE, indicates that we are performing a panel design and
#' the parameters used are specified in the panels and panel_overlap parameters.
#' @param number_panels The number of sample panels required.
#' @param verbose Boolean, set TRUE if you want to see output messaged to screen. Default is FALSE
#' i.e. no informational messages are displayed.
#'
#' @return Returns a list of the following variables:
#'
#' \itemize{
#'   \item \code{sample} This is a sample from the original shapefile that has had the
#'    appropriate panel id's add as a feature. The panel id values are determined by the
#'    panels and panel_overlap parameters.
#' }
#'
#' @keywords internal
PanelDesignAssignPanelids <- function(smp,
                                      panels,
                                      panel_overlap,
                                      panel_design,
                                      number_panels,
                                      verbose = FALSE){
  # if(panel_design) then assign panel id's to smp.
  # need to distinguish if panel_overlap is required.
  if(base::is.null(panel_overlap) & panel_design){
    if(verbose){
      base::message("spbal(PanelDesignAssignPanelids) Non-Overlapping panel design.")
    }
    tmp <- NULL
    # assign panel id's to sample points. smp$panel_id.
    for(i in 1:number_panels){
      tmp <- base::c(tmp, base::rep(i, panels[i]))
    }
    smp$panel_id <- tmp
  }
  # if(panel_overlap) is not null and panel_design is TRUE
  if(!base::is.null(panel_overlap) & panel_design){
    if(verbose){
      base::message("spbal(PanelDesignAssignPanelids) Overlapping panel design.")
    }
    # need to create the panel_id column
    # Initialize variables
    smp$panel_id <- 0
    panelid <- 1
    start_index <- 1

    for(i in 1:base::length(panels)){
      start_index <- start_index - panel_overlap[i]
      if(i == 1){
        for(j in 1:panels[i]){
          smp$panel_id[start_index] <- panelid
          start_index <- start_index + 1
        }
      } else {
        for(j in 1:panels[i]){
          smp$panel_id[start_index] <- base::list(base::c(smp$panel_id[start_index], panelid))
          start_index <- start_index + 1
        }
      }
      panelid <- panelid + 1
    }
  }
  result <- base::list(sample = smp)
  return(result)
}

