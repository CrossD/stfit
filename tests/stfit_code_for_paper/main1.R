################################
####### Code for paper #########
################################
library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(stfit)
library(RColorBrewer)

##################################################
############## Experiment setting ################
##################################################

#### Save all plots/tables to 'output' folder of the current working repository
if(!dir.exists("output"))
  dir.create("output")
#### Image theme setting
colthm = rasterTheme(panel.background=list(col="black"),region = brewer.pal(9, 'YlOrRd'))
colthm1 = RdBuTheme(panel.background=list(col="black"))
colthm1$regions$col = rev(colthm1$regions$col)
colthm2 = rasterTheme(panel.background=list(col="black"), region = "white")
#### Data setting
## The two datasets used in the paper is shipped with the package:
## - Site A: landsat2
## - Site B: landsat106
## Load data for the two sites
dfA = landsat2 %>% filter(year >= 2000)
dfB = landsat106 %>% filter(year >= 2000)
matA = as.matrix(dfA[,-c(1:2)])
matB = as.matrix(dfB[,-c(1:2)])
year = dfB$year
doy = dfB$doy

## Partial observed image index
pidx = c(68, 209, 352, 605, 624, 74, 156, 263, 273, 499, 184, 369, 508, 517, 565)
## Missing pattern matrix
pmat = matB[pidx,]
pmat[!is.na(pmat)] = 1
## Full observed image index
fidx1a = c(13, 101, 267, 432, 485)
fidx2a = c(21, 110, 192, 280, 493)
fidx3a = c(33, 121, 295, 458, 563) #609 
fidx4a = c(95, 128, 222, 261) 
fidxa = c(fidx1a, fidx2a, fidx3a, fidx4a)

fidx1b = c(145, 387, 481, 581, 587)
fidx2b = c(198, 276, 444, 493, 549)
fidx3b = c(82, 202, 293, 505, 557) #609 
fidx4b = c(132, 261, 265, 615, 657) 
fidxb = c(fidx1b, fidx2b, fidx3b, fidx4b)


#################################################
############## Images and tables ################
#################################################
#### Figure 1
pdf(paste0("output/fig_01.pdf"), width = 8, height = 5)
landsatVis(dfB[dfB$year == 2015, -c(1:2)], 
           names.attr = as.character(dfB$doy[dfB$year == 2015]),
           layout = c(9,5))
dev.off()

#### Fig 3
pdf("output/fig_03_left.pdf", width =6, height=5.5)
landsatVis(dfA[116, -c(1:2)], margin = FALSE)
dev.off()
pdf("output/fig_03_right.pdf", width =6, height=5.5)
landsatVis(dfB[198, -c(1:2)], margin = FALSE)
dev.off()

#### Fig 4
pdf("output/fig_04.pdf", width=6.9, height=5.5)
landsatVis(pmat, colthm = colthm2,
           names.attr = paste0("P", 1:15),
           layout = c(5, 3), colorkey=FALSE)
dev.off()

################ Illustration of the STFIT algorithm #################

#### 0. (Optional) Customize function for mean estimation ####
## Use 6 cores for parallel computating
registerDoParallel(3)
## Define a customized smoothing function using fourier basis
.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
## one dimentional smoothing regression function; given x and y return predicted y given x
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X[x.eval,] %*% lmfit$coefficient)
}
## Use 'customfun' function for mean estimation
## default is using 'smooth_spline' function
stfit::opts$set(temporal_mean_est = customfun)

## Create the artifically dataset 'mat' based on 'matB'
mat = matB
## Select the 198th image, which is fully observed as an example to illustrate the algorithm
fidx = 198
year[fidx] ## 2004
doy[fidx] ## 228
## Select the 8th missing pattern; midx: a logitcal vector of length 961, TRUE for missing pixel
midx = is.na(pmat[8,])
## missing pattern is applied to the fully observed image
mat[fidx, midx] = NA

## a list to save raster objects for plots later.
rst_lst=list()

## original image
rst_lst[[1]] = raster(matrix(matB[fidx, ], 31)) 
## image after applying the missing pattern
rst_lst[[2]] = raster(matrix(mat[fidx,], 31)) 

#### 1. Mean estimation ####
meanest = meanEst(doy, mat, doyeval = 1:365, clipRange = c(0,1800),
                  clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
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
## image for mean estimation
rst_lst[[3]] = raster(matrix(mat_mean_imp[fidx,], 31))

## calculate the residual matrix
rmat = mat - meanest$meanmat[unlist(lapply(doy, function(x,y) which(y == x), y = meanest$doyeval)),]
## residual image after removing mean
rst_lst[[4]] = raster(matrix(rmat[fidx,], 31))

#### 2. temporal effect estimation #### 
## teffarray is 16X365X961 array, where the first dimention is year
teffres = teffEst(year, doy, rmat, doyeval = meanest$doyeval, h.cov = 100, h.sigma2 = 300, var.est = TRUE)
teffarray = teffres$teff_array
teffvararray = teffres$teff_var_array
## matrix using mean + temporal effect for imputation; same dimention as mat
mat_teff_imp = mat_mean_imp
yearidx = unlist(lapply(year, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[1]])))
doyidx = unlist(lapply(doy, function(x,y) which(y == x), y = as.numeric(dimnames(teffarray)[[2]])))
for(i in 1:nrow(mat_teff_imp)){
  mat_teff_imp[i,] = mat_teff_imp[i,] + teffarray[yearidx[i], doyidx[i],]
}
## image for temporal effect estimation
rst_lst[[5]] = raster(matrix(teffarray[yearidx[fidx], doyidx[fidx],], 31))

#### 2.1 temporal effect illustration using pixel 157 #### 
## Visualization of temporal effect for pixel 157
pixelidx = 157
## Fig 8 Right
pdf("output/fig_08_right.pdf", width =6, height=5.5)
par(mar=c(4,4,0.1,0.1))
yeareval = as.numeric(dimnames(teffarray)[[1]])
doyeval = as.numeric(dimnames(teffarray)[[2]])
plot(doy, rmat[,pixelidx], col = gray(0.8), pch=19, ylab = "Residuals", xlab= "DOY")
abline(h=0, lty=2, col=gray(0.5))
ind = which(year==2004)
points(doy[ind], rmat[ind,pixelidx], col = 2, pch=19)
lines(doyeval, teffarray[5,,pixelidx], col = 2, lwd = 2)
dev.off()

## Visualization of mean + temporal effect for the target pixel
# pdf("output/fig_08_right1.pdf", width =6, height=5.5)
# par(mar=c(4,4,1,1))
# plot(doy, mat[,pixelidx], col = gray(0.8), pch=19, ylab = "y", xlab= "DOY")
# lines(doyeval, meanest$meanmat[,pixelidx], lwd = 2)
# ind = which(year==2004)
# points(doy[ind], mat[ind,pixelidx], col = 2, pch=19)
# lines(doyeval, meanest$meanmat[,pixelidx] + teffarray[5,,pixelidx], col = 2, lwd = 2, lty=2)
# dev.off()

weight.cov = weightVector(100)
yeareval = sort(unique(year))
t.grid = doyeval[unique(round(seq(1, length(doyeval), 
                                  length.out = 50)))]
resid = rmat[, pixelidx]
nnaidx = !is.na(resid)
## covariance function estimation
R0.hat = lc_cov_1d_est(year[nnaidx], doy[nnaidx], resid[nnaidx], 
                       weight.cov, t.grid)
## Fig 8 left
pdf("output/fig_08_left.pdf", width =6, height=5.5)
par(mar=c(0,0,0,0))
persp(R0.hat,theta=30, phi=30, expand=0.5, col='lightblue',
      xlab='DOY',ylab='DOY',zlab="Rs",ticktype='simple')
dev.off()


#### seffect estimation ####
## update residual matrix by removing temoral effect
for (i in 1:nrow(rmat)) {
  rmat[i, ] = rmat[i, ] - teffarray[yearidx[i], doyidx[i], ]
}
## Spatial effect estimation
seffres = seffEst(rmat, 31, 31, nnr = 30, h.cov = 2, h.sigma2 = 2, var.est = TRUE)
seffmat = seffres$seff_mat
seffvarmat = seffres$seff_var_mat
## matrix using mean + temporal effect + spatial effect for imputation; same dimention as mat
mat_seff_imp = mat_teff_imp + seffmat
## recover observed pixels
mat_seff_imp[fidx, !midx] = matB[fidx, !midx]

## residual image after removing temporal effect
rst_lst[[6]] = raster(matrix(rmat[fidx,], 31))

## residual median before adjusting for the temporal effect
median(values(rst_lst[[4]]), na.rm=TRUE)
## residual median after adjusting for the temporal effect
median(values(rst_lst[[6]]), na.rm=TRUE)

## estimated spatial effect image
rst_lst[[7]] = raster(matrix(seffmat[fidx,], 31))
## imputated image using mean + temporal effect + spatical effect
rst_lst[[8]] = raster(matrix(mat_seff_imp[fidx,], 31))
## residual matrix after further removing spatial effect
rmat = rmat - seffmat
## Final residual image
rst_lst[[9]] = raster(matrix(rmat[fidx,], 31))

## standard error
varest = teffvararray[yearidx[fidx], doyidx[fidx],] + seffvarmat[fidx,]
varest[!midx] = NA
rst_lst[[10]] = raster(sqrt(matrix(varest, 31)))

#### Plots

## Fig 5 Left
pdf("output/fig_05_left.pdf", width =6, height=5.5)
levelplot(rst_lst[[1]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
dev.off()
## Fig 5 Right
pdf("output/fig_05_right.pdf", width =6, height=5.5)
levelplot(rst_lst[[2]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
dev.off()
## Fig 6 Left
pdf("output/fig_06_left.pdf", width =6, height=5.5)
levelplot(rst_lst[[3]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
dev.off()
## Fig 6 Right
pdf("output/fig_06_right.pdf", width =6, height=5.5)
levelplot(rst_lst[[4]], par.settings = colthm1, at = seq(-245, 245, length.out = 20), margin = FALSE)
dev.off()
## Fig 7 Left
pdf("output/fig_07_left.pdf", width =6, height=5.5)
levelplot(rst_lst[[5]], par.settings = colthm1, at = seq(-245, 245, length.out = 20), margin = FALSE)
dev.off()
# pdf("output/fig_07_left1.pdf", width =6, height=5.5)
# levelplot(rst_lst[[3]] + rst_lst[[5]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
# dev.off()
## Fig 7 Right
pdf("output/fig_07_right.pdf", width =6, height=5.5)
levelplot(rst_lst[[6]], par.settings = colthm1, at = seq(-245, 245, length.out = 20), margin = FALSE)
dev.off()
## Fig 9 Left
pdf("output/fig_09_left.pdf", width =6, height=5.5)
levelplot(rst_lst[[7]], par.settings = colthm1, at = seq(-245, 245, length.out = 20), margin = FALSE)
dev.off()
## Fig 9 Right
pdf("output/fig_09_right.pdf", width =6, height=5.5)
levelplot(rst_lst[[9]], par.settings = colthm1, at = seq(-245, 245, length.out = 20), margin = FALSE)
dev.off()

# pdf("output/fig_10_left.pdf", width =6, height=5.5)
# levelplot(rst_lst[[3]] + rst_lst[[5]] + rst_lst[[7]], par.settings = colthm,
#           at = seq(200, 1300, length.out = 20), margin = FALSE)
# dev.off()
# # pdf("output/imputed.pdf", width =6, height=5.5)
# # levelplot(rst_lst[[8]], par.settings = colthm, at = seq(200, 1300, length.out = 20), margin = FALSE)
# # dev.off()
# pdf("output/fig_10_right.pdf", width =6, height=5.5)
# ## range(values(rst_lst[[10]]))
# levelplot(rst_lst[[10]], par.settings = colthm, margin = FALSE)
# dev.off()

################ Imputation  with stfit_landsat function ###############
pidx1 = pidx[c(6,8,14,15)]
pmat1 = pmat[c(6,8,14,15),]
#### fully observed image indexes from different seasons
fidx = fidxb[c(3,7,14,18)]
if(!dir.exists("output/tmp"))
  dir.create("output/tmp")

registerDoParallel(3)
rst_list1 = list()
for(i in 1:4){
  mat = matB
  ## apply missing patterns to fully observed images
  midx = is.na(pmat1[i,])
  mat[fidx[i], midx] = NA
  
  if(file.exists(paste0("output/tmp/res_P", pidx1[i], "_F", fidx[i], ".rds"))){
    res1 <- readRDS(paste0("output/tmp/res_P", pidx1[i], "_F", fidx[i], ".rds"))
  } else {
    res1 <- stfit_landsat(year, doy, mat, 31, 31, nnr=30,
                            use.intermediate.result = FALSE, intermediate.save = FALSE, var.est = TRUE)
    saveRDS(res1, paste0("output/tmp/res_P", pidx1[i], "_F", fidx[i], ".rds"))
  }
  # rst_list1[[(i-1)*3+1]] = raster(matrix(matB[fidx[i],], 31))
  # rst_list1[[(i-1)*3+2]] = raster(matrix(mat[fidx[i],], 31))
  # rst_list1[[(i-1)*3+3]] = raster(matrix(res1$imat[fidx[i],], 31))
  rst_list1[[(i-1)*4+1]] = raster(matrix(matB[fidx[i],], 31))
  rst_list1[[(i-1)*4+2]] = raster(matrix(mat[fidx[i],], 31))
  rst_list1[[(i-1)*4+3]] = raster(matrix(res1$imat[fidx[i],], 31))
  rst_list1[[(i-1)*4+4]] = raster(matrix(res1$sdmat[fidx[i],], 31))
}

s = stack(rst_list1)
pdf("output/fig_10.pdf", width =6, height=5.5)
# levelplot(s,par.settings = colthm, index.cond=list(c(1,4,7,10,2,5,8,11,3,6,9,12)),
#           names.attr = c(rbind(paste0("F", c(3,7,14,18)), paste0("F", c(3,7,14,18), "P", c(6,8,14,15)),
#                                paste0("F", c(3,7,14,18), "P", c(6,8,14,15), " impu"))))
# levelplot(s, par.settings = colthm, index.cond=list(c(seq(1, 16, 4), seq(2, 16, 4), seq(3, 16, 4), seq(4, 16, 4))),
#           names.attr = c(rbind(paste0("F", c(3,7,14,18)),
#                                paste0("F", c(3,7,14,18), "P", c(6,8,14,15)),
#                                paste0("F", c(3,7,14,18), "P", c(6,8,14,15), " impu"),
#                                paste0("F", c(3,7,14,18), "P", c(6,8,14,15), " sd"))),
#           layout = c(4,4), zscaleLog = TRUE)
levelplot(s, par.settings = rasterTheme(panel.background=list(col="black"),
                                        region = brewer.pal(9, 'YlOrRd')[1:7]),
          index.cond=list(c(seq(1, 16, 4), seq(2, 16, 4), seq(3, 16, 4), seq(4, 16, 4))),
          names.attr = c(rbind(paste0("F", c(3,7,14,18)),
                               paste0("F", c(3,7,14,18), "P", c(6,8,14,15)),
                               paste0("F", c(3,7,14,18), "P", c(6,8,14,15), " impu"),
                               paste0("F", c(3,7,14,18), "P", c(6,8,14,15), " sd"))),
          layout = c(4,4), zscaleLog = TRUE)

dev.off()

# levelplot(stack(rst_list1[seq(1, 16, 4)]), par.settings = colthm, layout=c(4,1),
#           names.attr = paste0("F", c(3,7,14,18)))
# levelplot(stack(rst_list1[seq(2, 16, 4)]), par.settings = colthm, layout=c(4,1),
#           names.attr = paste0("F", c(3,7,14,18), "P", c(6,8,14,15)))
# pdf("output/fig_12.pdf", width =6, height=5.5)
# levelplot(stack(rst_list1[seq(4, 16, 4)]), par.settings = colthm, layout=c(4,1),
#           names.attr = paste0("F", c(3,7,14,18), "P", c(6,8,14,15), " sd"))
# dev.off()
# levelplot(stack(rst_list1[c(seq(2, 16, 4), seq(4, 16, 4))]), par.settings = colthm,
#           layout=c(4,2))

sapply(seq(4, 16, 4), function(i) sd(values(rst_list1[[i]])))

################ Supplementary #################

## Fig S1
pdf("output/fig_S1.pdf", width=6.9, height=5.5)
landsatVis(matA[fidxa,], 
           names.attr = paste0("F", 1:19), layout = c(5,4))
dev.off()
## Fig S2 
pdf("output/fig_S2.pdf", width=6.9, height=5.5)
landsatVis(matB[fidxb,], 
           names.attr = paste0("F", 1:20), layout = c(5,4))
dev.off()
## Table S1 (a)
print(xtable::xtable(
  data.frame(ID = paste0("F", 1:19),
             year = dfA$year[fidxa],
             DOY = dfA$doy[fidxa])
), include.rownames=FALSE, file = "output/tab_S1a.txt")

## Table S1 (b)
print(xtable::xtable(
  data.frame(ID = paste0("F", 1:20),
             year = dfB$year[fidxb],
             DOY = dfB$doy[fidxb])
), include.rownames=FALSE, file = "output/tab_S1b.txt")

## Table S2
print(xtable::xtable(
  data.frame(ID = paste0("P", 1:15),
             pct=apply(pmat, 1,
                       function(x) round(sum(is.na(x))/length(x), 2)))), 
  include.rownames=FALSE, file = "output/tab_S2.txt")



