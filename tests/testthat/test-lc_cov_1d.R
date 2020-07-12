test_that("lc_cov_1d function check", {
  # skip("skip this test")
  skip_on_cran()
  lc_cov_1d_r <- function(ids, time, resid, W, t1,t2){
    W_size = length(W)
    sumEEKK = 0.0
    sumKK = 0.0
    N = length(ids)
    time_min = min(time);
    time_max = max(time);
    
    k1_start = max(t1 - W_size%/%2+1, time_min);
    k2_start = max(t2 - W_size%/%2+1, time_min);
    k1_stop = min(t1 + W_size%/%2+1, time_max);
    k2_stop = min(t2 + W_size%/%2+1, time_max);
    
    for (i in 1:N){
      # if(i %in% c(32, 48, 63, 87, 158, 184))
      #   browser()
      if(time[i] >= k1_start & time[i] < k1_stop){
        for(j in 1:N){
          # if(j %in% c(32, 48, 63, 87, 158, 184))
          #   browser()
          if(i == j)
            next
          if(ids[i] == ids[j]){
            if(time[j] >= k2_start & time[j] < k2_stop){
              sumEEKK = sumEEKK + resid[i]*resid[j]*W[time[i] - t1 + W_size%/%2+1]*W[time[j] - t2 + W_size%/%2+1];
              sumKK = sumKK +  W[time[i] - t1 + W_size%/%2+1]*W[time[j] - t2 + W_size%/%2+1];
            }
          }
        }
      }
    }
    # browser()
    if(sumKK == 0.0){
      return(NA);
    } else{
      return(sumEEKK/sumKK);
    }
  }
  
  ## test
  tstdat = readRDS(system.file("testdata", "cov_1d_test_data.rds", package = "stfit"))
  with(tstdat, {
    h = 100
    W = epan(seq(-h, h, by = 1)/h)
    expect_equal(stfit:::lc_cov_1d(ids, x, resid, W, 1,1), lc_cov_1d_r(ids, x, resid, W,1,1))
  })
})
