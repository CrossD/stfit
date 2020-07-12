library(feather)
library(dplyr)
library(geoR)
library(raster)
library(rasterVis)

df = landsat2 %>% filter(year >= 2000)
year = df$year
doy = df$doy
mat0 = as.matrix(df[,-c(1:2)])

#######################################
##### Simulation study for spring #####
#######################################
#### partial missing image indexes with different missing percentage
## Partial observed image index
pidx = c(68, 209, 352, 605, 624, 74, 156, 263, 273, 499, 184, 369, 508, 517, 565)
pmat = as.matrix(landsat106[landsat106$year >= 2000,-c(1:2)])[pidx,]
pmat[!is.na(pmat)] = 1
#### fully observed image indexes from different seasons
fidx1 = c(13, 101, 267, 432, 485)
fidx2 = c(21, 110, 192, 280, 493)
fidx3 = c(33, 121, 295, 458, 563) #609 
fidx4 = c(95, 128, 222, 261) 
fidx = c(fidx1, fidx2, fidx3, fidx4)
fmat = mat0[fidx, ]
if(!dir.exists("output"))
  dir.create("output")

###### Detrend images with the same mean estimation procedure as STFIT first #####
###### define function used for mean estimation ######
.X = fda::eval.basis(1:365, fda::create.fourier.basis(rangeval=c(0,365), nbasis=11))
customfun <- function(x, y, x.eval=1:365, minimum.num.obs = 10){
  nonna.idx = !is.na(y)
  if(sum(nonna.idx) < minimum.num.obs)
    return(rep(NA, 365))
  ## lmfit = lm.fit(.X[unlist(lapply(x, function(x) which(x == x.eval))),], y[nonna.idx])
  lmfit = lm.fit(.X[x[nonna.idx],], y[nonna.idx])
  return(.X[x.eval,] %*% lmfit$coefficient)
}
stfit::opts$set(temporal_mean_est = customfun)
registerDoParallel(3)

RMSEmat = matrix(NA, nrow(pmat), length(fidx))
NMSEmat = matrix(NA, nrow(pmat), length(fidx))
AREmat = matrix(NA, nrow(pmat), length(fidx))
CORmat = matrix(NA, nrow(pmat), length(fidx))
for(i in 1:nrow(pmat)){
  for(j in 1:length(fidx)){
    fmatj = fmat[j,]
    ## apply missing patterns to fully observed images
    missing.idx = is.na(pmat[i,])
    mat[fidx[j], missing.idx] = NA
    fmatj[missing.idx] = NA
    
    ## detrend by using mean estimation
    meanest = meanEst(doy, mat, doyeval = 1:365, outlier.tol = 0.2, clipRange = c(0, 1800),
                      clipMethod = "nnr", img.nrow = 31, img.ncol = 31)
    ## remove outlier pixels
    if(fidx[j] %in% meanest$outlier$outidx)
      fmatj[meanest$outlier$outlst[[which(fidx[j] == meanest$outlier$outidx)]]] = NA
    rmatj = fmatj - meanest$meanmat[which(doy[fidx[j]] == meanest$doyeval),]
    
    gdata = as.geodata(cbind(expand.grid(seq(1,31), seq(1, 31)), rmatj))
    vario = variog(gdata)
    wls = variofit(vario, ini = c(150000, 35))
    loci <- expand.grid(seq(0,1,l=31), seq(0,1,l=31))
    # predicting by ordinary kriging
    kc <- krige.conv(gdata, loc=expand.grid(seq(1,31), seq(1, 31)),
                     krige=krige.control(cov.pars= wls$cov.pars))
    imat = meanest$meanmat[which(doy[fidx[j]] == meanest$doyeval),] + kc$predict
    # image(kc, main="kriging estimates")
    # image(matrix(fmatj, 31))
    #### proposed method
    RMSEmat[i, j] = RMSE(fmat[j, missing.idx], imat[missing.idx])
    NMSEmat[i, j] = NMSE(fmat[j, missing.idx], imat[missing.idx])
    AREmat[i, j] = ARE(fmat[j, missing.idx], imat[missing.idx])
    CORmat[i, j] = cor(fmat[j, missing.idx], imat[missing.idx])
  }
}

saveRDS(list(RMSEmat, NMSEmat, AREmat, CORmat), "./output/res.rds")

