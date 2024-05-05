# Validate Halton Frame functions, features and parameter validation.
# test-spbal-HF-1.R

testthat::test_that("1. Verify panels= and panel_overlap= parm length checking.", {
  n_panels <- c(20, 20, 20, 100)
  n_panel_overlap <- c(0, 4, 5)
  shp.cant <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  bb <- spbal::BoundingBox(shp.cant)
  expect_error(spbal::HaltonFrame(shapefile = shp.cant,
                                  panels = n_panels,
                                  panel_overlap = n_panel_overlap,
                                  boundingbox = bb), "spbal(ValidatePanelDesign) length of panels [4] must match length of panel_overlap [3].", fixed=TRUE)
})

testthat::test_that("2. Verify panels= and panel_overlap= parm length checking.", {
  n_panels <- c(20, 20, 20)
  n_panel_overlap <- c(0, 4, 5, 6)
  shp.cant <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  bb <- spbal::BoundingBox(shp.cant)
  expect_error(spbal::HaltonFrame(shapefile = shp.cant,
                                  panels = n_panels,
                                  panel_overlap = n_panel_overlap,
                                  boundingbox = bb), "spbal(ValidatePanelDesign) length of panels [3] must match length of panel_overlap [4].", fixed=TRUE)
})

testthat::test_that("3. Ensure first HaltonFrame point has SiteID == 1.", {
  set.seed(511)
  # read sample shapefile from sf package.
  sf_object <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  # Vertically aligned master sample bounding box.
  bb <- spbal::BoundingBox(shapefile = sf_object)
  # generate 100 samples
  n_samples <- 100
  # suppress "attribute variables are assumed to be spatially constant throughout all geometries"
  suppressWarnings(result <- spbal::HaltonFrame(shapefile = sf_object,
                                                N = n_samples,
                                                boundingbox = bb,
                                                verbose = FALSE) )
  HF100 <- result$hf.pts.shp
  expect_equal(HF100[1,]$ID, 1)
})

testthat::test_that("4. Ensure first HaltonFrame point has spbalSeqID == 1.", {
  # read sample shapefile from sf package.
  sf_object <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  # Vertically aligned master sample bounding box.
  bb <- spbal::BoundingBox(shapefile = sf_object)
  # generate 100 samples
  n_samples <- 100
  # suppress "attribute variables are assumed to be spatially constant throughout all geometries"
  suppressWarnings(result <- spbal::HaltonFrame(shapefile = sf_object,
                                                N = n_samples,
                                                boundingbox = bb,
                                                verbose = FALSE) )
  HF100 <- result$hf.pts.shp
  expect_equal(HF100[1,]$spbalSeqID, 1)
})

