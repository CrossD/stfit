#include <Rcpp.h>
#include <stdlib.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
double vecmin(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}
// [[Rcpp::export]]
double vecmax(NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}
double lc_cov_1d_(const NumericVector &ids, const NumericVector &time, const NumericVector &resid, 
                  const NumericVector &W, int t1, int t2){
  // sparse local constant covariance estimation for points i and j
  // ids: vector of ids
  // time: vector of observed time points
  // resid: vector of residuals
  // W: weight vector, related to the bandwidth selection
  // t1: first time point value; t2: second time point value
  int W_size = W.size();
  double sumEEKK = 0.0, sumKK = 0.0;
  int N = ids.size();
  int time_min = (int)vecmin(time);
  int time_max = (int)vecmax(time);

  /* the starts */
  int k1_start = std::max(t1 - W_size/2, time_min);
  int k2_start = std::max(t2 - W_size/2, time_min);
  
  /* the stops */
  int k1_stop = std::min(t1 + W_size/2 + 1, time_max);
  int k2_stop = std::min(t2 + W_size/2 + 1, time_max);

  for(int i = 0; i < N; i++){
    if(time[i] >= k1_start & time[i] < k1_stop){
      for(int j = 0; j < N; j++){
        if(i == j)
          continue;
        if(ids[i] == ids[j]){
          if(time[j] >= k2_start & time[j] < k2_stop){
            sumEEKK += resid[i]*resid[j]*W[time[i] - t1 + W_size/2]*W[time[j] - t2 + W_size/2];
            sumKK += W[time[i] - t1 + W_size/2]*W[time[j] - t2 + W_size/2];
          }
        }
      }
    }
  }
  if(sumKK == 0.0){ // no points within the bandwidth
    Rcpp::Rcout << "sumKK is 0" << std::endl;
    return NA_REAL;
  } else{
    return sumEEKK/sumKK;
  }
}

// [[Rcpp::export]]
double lc_cov_1d(const NumericVector &ids, const NumericVector &time, const NumericVector &resid, 
                 const NumericVector &W, int t1, int t2){
  return lc_cov_1d_(ids, time, resid, W, t1, t2);
}


// [[Rcpp::export]]
NumericMatrix lc_cov_1d_est(const NumericVector &ids, const NumericVector &time, const NumericVector &resid, 
                  const NumericVector &W, const NumericVector &tt){
  // X.nrow() == nRow * nCol
  // each row of X is a row stacked image
  int tt_size = tt.size();
  NumericMatrix out(tt_size, tt_size);

  for(int i = 0; i < tt_size; i++){
    // if(i == 0)
    //   Rcpp::Rcout << "tt[i]: " << tt[i] << std::endl;
    for(int j = 0; j <= i; j++){
       out(i,j) = lc_cov_1d_(ids, time, resid, W, tt[i], tt[j]);
      if(j < i)
        out(j,i) = out(i,j);
      // if(NumericVector::is_na(out(i,j)))
      //   Rcpp::Rcout << "i: " << i << "; j: " << j<< std::endl;
    }
  }
  return(out);
}