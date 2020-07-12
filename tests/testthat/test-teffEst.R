library(testthat)
library(stfit)
library(dplyr)
library(doParallel)
registerDoParallel(3)
test_that("teffEst backward compatibility check", {
  # skip("skip this test")
  skip_on_cran()
  dfB = landsat106 %>% filter(year >= 2000)
  year = dfB$year
  doy = dfB$doy
  rmat = readRDS(system.file("testdata", "rmat1_B.rds", package = "stfit"))
  teffarray = teffEst(year, doy, rmat, doyeval = 1:365, h.cov = 100, h.sigma2 = 300)$teff_array
  ## saveRDS(teffarray, "../../inst/testdata/teffarray_B.rds")
  expect_identical(teffarray, readRDS(system.file("testdata", "teffarray_B.rds", package = "stfit")))
})
