library(feather)
library(dplyr)
library(doParallel)
library(Matrix)
library(raster)
library(rasterVis)
library(fda)
library(stfit)
library(abind)

df = landsat106 %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])

#######################################
##### Simulation study for spring #####
#######################################
#### partial missing image indexes with different missing percentage
##### selected partially observed images indexes
pidx = c(68, 209, 352, 605, 624, 74, 156, 263, 273, 499, 184, 369, 508, 517, 565)
pmat = as.matrix(landsat106[landsat106$year >= 2000,-c(1:2)])[pidx,]
pmat[!is.na(pmat)] = 1
#### fully observed image indexes from different seasons
fidx1 = c(145, 387, 481, 581, 587)
fidx2 = c(198, 276, 444, 493, 549)
fidx3 = c(82, 202, 293, 505, 557) #609 
fidx4 = c(132, 261, 265, 615, 657) 
fidx = fidx1
fmat = mat0[fidx, ]
if(!dir.exists("stfit_spring"))
  dir.create("stfit_spring")

###### define function used for mean estimation ######
.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X[x.eval,] %*% lmfit$coefficient)
}
stfit::opts$set(temporal_mean_est = customfun)

## collapse the partial and full index for parallel
## matrix of MxN, column stacking
N = nrow(fmat)
M = nrow(pmat)
registerDoParallel(3)
res = foreach(n = 1:(M*N)) %dopar% {
  i = (n - 1) %% M + 1 ## ROW INDEX
  j = (n - 1) %/% M + 1 ## COLUMN INDEX
  mat = mat0
  ## apply missing patterns to fully observed images
  missing.idx = is.na(pmat[i,])
  mat[fidx[j], missing.idx] = NA
  
  if(file.exists(paste0("./stfit_spring/stfit_spring_P", pidx[i], "_F", fidx[j], ".rds"))){
      res1 <- readRDS(paste0("./stfit_spring/stfit_spring_P", pidx[i], "_F", fidx[j], ".rds"))
  } else {
      res1 <- gapfill_landsat(year, doy, mat, 31, 31, nnr=30,
                              use.intermediate.result = FALSE, intermediate.save = FALSE)
      saveRDS(res1, paste0("./stfit_spring/stfit_spring_P", pidx[i], "_F", fidx[j], ".rds"))
  }
  imat = res1$imat[fidx[j],]
  c(RMSE(fmat[j, missing.idx], imat[missing.idx]),
    NMSE(fmat[j, missing.idx], imat[missing.idx]),
    ARE(fmat[j, missing.idx], imat[missing.idx]),
    cor(fmat[j, missing.idx], imat[missing.idx]))
}
saveRDS(res, "./stfit_spring/res.rds")

