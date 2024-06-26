% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getSample.R
\name{getSample}
\alias{getSample}
\title{Extract a sample of a specified size from a master sample.}
\usage{
getSample(shapefile, n, randomStart = FALSE, strata = NULL, stratum = NULL)
}
\arguments{
\item{shapefile}{A MULTIPOINT or POINT object from where to take the sample.}

\item{n}{The number of sample points to return.}

\item{randomStart}{Whether a spatially balanced sample will be randomly drawn from
the frame or not. Default value is FALSE.}

\item{strata}{to be added}

\item{stratum}{The name of a column in the dataframe attached to shapefile that defines
the strata of interest.}
}
\value{
A list containing the following variable:

\itemize{
\item \code{sample} The sample from the shapefile POINTS.
}
}
\description{
A description of this useful function.
}
\examples{
# Draw a spatially balanced sample of n = 25 from a Halton Frame over Gates --

# Use the North Carolina shapefile supplied in the sf R package.
shp_file <- sf::st_read(system.file("shape/nc.shp", package="sf"))
shp_gates <- shp_file[shp_file$NAME == "Gates",]

# Vertically aligned master sample bounding box.
bb <- spbal::BoundingBox(shapefile = shp_gates)

set.seed(511)
result7 <- spbal::HaltonFrame(shapefile = shp_gates,
                              J = c(6, 4),
                              boundingbox = bb)
Frame <- result7$hf.pts.shp

# Get the first 25 sites from a B = (2^6) * (3^4) Halton Frame (62,208 grid
# points covering Gates).
n_samples <- 25
FrameSample <-getSample(shapefile = Frame,
                        n = n_samples)
FrameSample <- FrameSample$sample
FrameSample

}
\author{
Phil Davies.
}
