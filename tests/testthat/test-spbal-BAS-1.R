# Validate BAS functions, features and parameter validation.
# test-spbal-BAS-1.R

testthat::test_that("1. Verify panels= and panel_overlap= parm length checking.", {
  n_panels <- c(20, 20, 20, 100)
  n_panel_overlap <- c(0, 4, 5)
  shp.nc <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  bb <- spbal::BoundingBox(shp.nc)
  expect_error(spbal::BAS(shapefile = shp.nc,
                          panels = n_panels,
                          panel_overlap = n_panel_overlap,
                          boundingbox = bb), "spbal(ValidatePanelDesign) length of panels [4] must match length of panel_overlap [3].", fixed=TRUE)
})

testthat::test_that("2. Verify panels= and panel_overlap= parm length checking.", {
  n_panels <- c(20, 20, 20)
  n_panel_overlap <- c(0, 4, 5, 6)
  shp.nc <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  bb <- spbal::BoundingBox(shp.nc)
  expect_error(spbal::BAS(shapefile = shp.nc,
                          panels = n_panels,
                          panel_overlap = n_panel_overlap,
                          boundingbox = bb), "spbal(ValidatePanelDesign) length of panels [3] must match length of panel_overlap [4].", fixed=TRUE)
})

testthat::test_that("3. Ensure first BAS point has SiteID == 1.", {
  # read sample shapefile from sf package.
  sf_object <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  # Vertically aligned master sample bounding box.
  bb <- spbal::BoundingBox(shapefile = sf_object)
  # generate 100 samples
  n_samples <- 100
  result <- spbal::BAS(shapefile = sf_object,
                       n = n_samples,
                       boundingbox = bb,
                       verbose = FALSE)
  BAS100 <- result$sample
  expect_equal(BAS100[1,]$SiteID, 1)
})

testthat::test_that("4. Ensure first BAS point has spbalSeqID == 1.", {
  # read sample shapefile from sf package.
  sf_object <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  # Vertically aligned master sample bounding box.
  bb <- spbal::BoundingBox(shapefile = sf_object)
  # generate 100 samples
  n_samples <- 100
  result <- spbal::BAS(shapefile = sf_object,
                       n = n_samples,
                       boundingbox = bb,
                       verbose = FALSE)
  BAS100 <- result$sample
  expect_equal(BAS100[1,]$spbalSeqID, 1)
})

testthat::test_that("5. Verify messaging when no shapefile specified.", {
  n_panels <- c(20, 20, 20, 100)
  n_panel_overlap <- c(0, 4, 5)
  shp.nc <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  bb <- spbal::BoundingBox(shp.nc)
  expect_error(spbal::BAS(panels = n_panels,
                          panel_overlap = n_panel_overlap,
                          boundingbox = bb), "spbal(BAS) The shapefile parameter must be used. Please specify a shapefile.", fixed=TRUE)
})



