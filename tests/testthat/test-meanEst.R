library(testthat)
library(stfit)
library(dplyr)
library(doParallel)
registerDoParallel(3)
test_that("meanEst backward compatibility check", {
  # skip("skip this test")
  skip_on_cran()
  dfB = landsat106 %>% filter(year >= 2000)
  mat = as.matrix(dfB[,-c(1:2)])
  year = dfB$year
  doy = dfB$doy
  
  meanest = meanEst(doy, mat, doyeval = 1:365, clipRange = c(0,1800),
                    clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
  expect_identical(meanest, readRDS(system.file("testdata", "meanest_B.rds", package = "stfit")))
})
