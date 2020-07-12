library(testthat)
library(stfit)
library(dplyr)
library(doParallel)
library(Matrix)
registerDoParallel(3)
test_that("seffEst backward compatibility check", {
  # skip("skip this test")
  skip_on_cran()
  rmat = readRDS(system.file("testdata", "rmat2_B.rds", package = "stfit"))
  seff_mat = seffEst(rmat, 31, 31, nnr = 30, h.cov = 2, h.sigma2 = 2)$seff_mat
  tmp = readRDS(system.file("testdata", "seffest_B.rds", package = "stfit"))
  expect_equal(sum(!tmp$seffmat[tmp$idx$idx.partialmissing,] == seff_mat[tmp$idx$idx.partialmissing,]), 0)
})
