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
  n_panels <- c(20, 20, 20)
  n_panel_overlap <- c(0, 4, 5)
  shp.nc <- sf::st_read(base::system.file("shape/nc.shp", package="sf"))
  bb <- spbal::BoundingBox(shp.nc)
  expect_error(spbal::BAS(panels = n_panels,
                          panel_overlap = n_panel_overlap,
                          boundingbox = bb), "spbal(BAS) The shapefile parameter must be used. Please specify a shapefile.", fixed=TRUE)
})

testthat::test_that("6. Verify BAS is drawing the exact same points from previous versions.", {
  ## Values drawn manually from original checked version main-branch of spbal on JavaLake2023 (5-13-2024).
  vals <- data.frame(X = c(641615.77490234375, 454764.05615234375, 604245.43115234375, 529504.74365234375, 536511.68310546875),
                     Y = c(906139.560585276805795729160308837890625, 927838.918609968037344515323638916015625,
                           341956.251943301293067634105682373046875, 537250.474165523657575249671936035156250,
                           732544.696387745789252221584320068359375))


  ## Draw the same BAS Sample and compare:
  bb.df <- c("xmin" = 85148, "ymin" = 33745, "xmax" = 1280999, "ymax" = 1351981)
  bb <- sf::st_as_sfc(sf::st_bbox(bb.df))
  msproj <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
  sf::st_crs(bb) <- msproj
  shp.df <- data.frame(id = 1:3, x = c(408891.353328072465956211090087890625, 965625.853148331865668296813964843750, 543869.399493260309100151062011718750),
                       y = c(933567.004960335791110992431640625000, 991719.841984464786946773529052734375, 85858.92871646676212549209594726562500))
  shp <- sf::st_as_sf(shp.df, coords = c("x", "y"), crs = msproj)
  shp <- sf::st_cast(sf::st_combine(shp), "POLYGON")

  ## Run BAS and get first 5 points:
  smp <- spbal::BAS(shapefile = shp, n = 5, boundingbox = bb, seeds = c(1257,3557))
  coords <- as.data.frame(sf::st_coordinates(smp$sample))

  ## Throw an error if any updates in the BAS code change these points.
  expect_equal(coords, vals)
})



