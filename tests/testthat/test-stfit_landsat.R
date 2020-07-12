library(testthat)
library(stfit)
library(dplyr)
library(doParallel)
library(Matrix)
registerDoParallel(3)
test_that("stfit_landsat backward compatibility check", {
  # skip("skip this test")
  skip_on_cran()
  dfB = landsat106 %>% filter(year >= 2000)
  mat = as.matrix(dfB[,-c(1:2)])
  year = dfB$year
  doy = dfB$doy
  res <- stfit_landsat(year, doy, mat, 31, 31, nnr=30, var.est = TRUE)
  tmp = readRDS(system.file("testdata", "seffest_B.rds", package = "stfit"))
  expect_equal(res$imat, tmp$imat)
  expect_equal(res$sdmat, tmp$sdmat)
})
