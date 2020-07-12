## Code used to generate test data
library(dplyr)
dfB = landsat106 %>% filter(year >= 2000)
mat = as.matrix(dfB[,-c(1:2)])
year = dfB$year
doy = dfB$doy
registerDoParallel(3)

#########################################
## test data for lc_cov_1d* functions ###
#########################################
y = matB[,3]
nnaidx = !is.na(y)
ids = df$year[nnaidx]
x = df$doy[nnaidx]
y = y[nnaidx]
uni.ids = unique(ids)

## mean estimation with smoothing spline
yfit = smooth_spline(x, y)
resid = y - yfit
.outlier <- function(y){
  whisker = boxplot(y, plot = FALSE)$stats[c(1, 5)]
  which(y < whisker[1] | y > whisker[2])
}
outid = .outlier(resid)

## remove outliers
ids = ids[-outid]
x = x[-outid]
y = y[-outid]

## individual curves
# plot(x, y, ylim = range(y), type = "n")
# ##for(i in 1:length(uni.ids)){
# for(i in 1:length(uni.ids)){
#   idx = ids == uni.ids[i]
#   lines(x[idx], y[idx])
# }
## redo mean estimation
yfit = smooth_spline(x, y)
# lines(sort(x), yfit[order(x)], lwd = 4)

## residual plots
resid = y - yfit
# plot(x, resid, ylim = range(resid), type = "n")
# ##for(i in 1:length(uni.ids)){
# for(i in 1:length(uni.ids)){
#   idx = ids == uni.ids[i]
#   lines(x[idx], resid[idx])
# }

cov_1d_test_data = list(ids=ids, x=x, resid=resid)
saveRDS(cov_1d_test_data, file = "inst/testdata/cov_1d_test_data.rds")


######################################################
## test data for meanEst/teffEst/seffEst functions ###
######################################################
meanest = meanEst(doy, mat, doyeval = 1:365, clipRange = c(0,1800),
                  clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
saveRDS(meanest, file = "inst/testdata/meanest_B.rds")

## remove outlier pixels
for(i in 1:length(meanest$outlier$outidx)){
  mat[meanest$outlier$outidx[i], meanest$outlier$outlst[[i]]] = NA
}
## remove outlier images
outlier.img.idx = meanest$idx$idx.outlier
for(i in outlier.img.idx){
  mat[outlier.img.idx,] = NA
}
## mean estimation after removing outliers
meanest = meanEst(doy, mat, doyeval = 1:365, clipRange = c(0,1800),
                  clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
## matrix using mean estimation for imputation; same dimention as mat
mat_mean_imp = meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]

## calculate the residual matrix
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
saveRDS(rmat, file = "inst/testdata/rmat1_B.rds")

#### 2. temporal effect estimation #### 
## teffarray is 16X365X961 array, where the first dimention is year
teffarray = teffEst(year, doy, rmat, doyeval = meanest$doyeval, h.cov = 100, h.sigma2 = 300)
saveRDS(teffarray, file = "inst/testdata/teffarray_B.rds")

#### seffect estimation ####
## update residual matrix by removing temoral effect
for (i in 1:nrow(rmat)) {
  rmat[i, ] = rmat[i, ] - teffarray[yearidx[i], doyidx[i], ]
}
saveRDS(rmat, file = "inst/testdata/rmat2_B.rds")

## Spatial effect estimation
seffest = seffEst(rmat, 31, 31, nnr = 30, h.cov = 2, h.sigma2 = 2)
saveRDS(seffest, file = "inst/testdata/seffest_B.rds")

###########################################
## test data for stfit_landsat function ###
###########################################
res <- stfit_landsat(year, doy, mat, 31, 31, nnr=30,
use.intermediate.result = FALSE, intermediate.save = FALSE, var.est = TRUE)
saveRDS(res, file = "inst/testdata/stfit_landsat_B.rds")



