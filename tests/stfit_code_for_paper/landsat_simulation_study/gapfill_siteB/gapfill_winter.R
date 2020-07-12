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
##### Simulation study for winter #####
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
fidx = fidx4
fmat = mat0[fidx, ]
if(!dir.exists("gapfill_winter"))
  dir.create("gapfill_winter")


###### variables used for Gapfill package ######
doybin = findInterval(doy, seq(1,365, by=8))
yearuni = sort(unique(year))
doybinuni = sort(unique(doybin))

## collapse the partial and full index for parallel
## matrix of MxN, column stacking
N = nrow(fmat)
M = nrow(pmat)
registerDoParallel(10)
res = foreach(n = 1:(M*N)) %dopar% {
  i = (n - 1) %% M + 1 ## ROW INDEX
  j = (n - 1) %/% M + 1 ## COLUMN INDEX
  mat = mat0
  ## apply missing patterns to fully observed images
  missing.idx = is.na(pmat[i,])
  mat[fidx[j], missing.idx] = NA
  
  datarray = array(NA, dim = c(31, 31, 46, 16), dimnames = list(1:31, 1:31, doybinuni, yearuni))
  for(ii in 1:16){
    for(jj in 1:46){
      idx = year == yearuni[ii] & doybin == doybinuni[jj]
      if(sum(idx) == 1)
        datarray[,,jj,ii] = matrix(mat[year == yearuni[ii] & doybin == doybinuni[jj],], 31) else
          if(sum(idx) > 1)
            warning("Multiple matches.")
    }
  }
  yidx = which(year[fidx[j]] == yearuni)
  didx = which(findInterval(doy[fidx[j]], seq(1,365, by=8)) == doybinuni)
  didxinterval = max(1,didx-6):min(46, didx + 6)
  yidxinterval = max(1, yidx - 4):min(16, yidx + 4)
  tmpmat = datarray[,,didxinterval, yidxinterval]
  
  if(file.exists(paste0("./gapfill_winter/gapfill_winter_P", pidx[i], "_F", fidx[j], ".rds"))){
    res1 <- readRDS(paste0("./gapfill_winter/gapfill_winter_P", pidx[i], "_F", fidx[j], ".rds"))
  } else {
    res1 = gapfill::Gapfill(tmpmat, clipRange = c(0, 1800), dopar = TRUE)
    saveRDS(res1, paste0("./gapfill_winter/gapfill_winter_P", pidx[i], "_F", fidx[j], ".rds"))
  }
  imat = c(res1$fill[,,which(didx == didxinterval), which(yidx == yidxinterval)])
  c(RMSE(fmat[j, missing.idx], imat[missing.idx]),
    NMSE(fmat[j, missing.idx], imat[missing.idx]),
    ARE(fmat[j, missing.idx], imat[missing.idx]),
    cor(fmat[j, missing.idx], imat[missing.idx]))
}
saveRDS(res, "./gapfill_winter/res.rds")
