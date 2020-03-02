/*!
  @file cumres.cpp
  @author Klaus K. Holst
  @copyright 2020, Klaus Kähler Holst

  @brief Generic function for calculating cumulative residuals

  Test statistics 

*/

#include "cumres.hpp"

namespace target {

  cumres::cumres(const arma::vec &r, const arma::mat &dr, const arma::mat &ic) : r(r), dr(dr), ic(ic) {
#ifndef ARMA_R
    arma::arma_rng::set_seed_random();
#endif
    n = r.n_elem;
    arma::vec inp(n);
    for (unsigned i=0; i<n; i++) inp(i) = i;
    this->ord = arma::conv_to<arma::uvec>::from(inp);    
    this->reorder(inp);
  }
 
  void cumres::reorder(const arma::vec &inp) {    
    arma::uvec revord = arma::stable_sort_index(ord); // back to original order of input data
    ord = arma::stable_sort_index(inp); // new order
    arma::vec tt = t;
    t = inp.elem(ord);
    revord = revord.elem(ord);
    r = r.elem(revord);
    dr = dr.rows(revord);
    eta = arma::cumsum(dr, 0); // cumulative sum of each column 
    ic = ic.rows(revord);
  }

  arma::vec cumres::rnorm() {
#ifdef ARMA_R
    Rcpp::RNGScope scope;
    return Rcpp::as<arma::vec>(Rcpp::rnorm(n));    
#else
    return arma::randn<arma::vec>(n);
#endif    
  }
  
  arma::vec cumres::obs() {
    return arma::cumsum(r)/std::sqrt((double)n);
  }

  // Sample single process 
  arma::vec cumres::sample(arma::uvec idx) { 
    unsigned N = n;
    if (!idx.is_empty()) {
      N = idx.n_elem;
    }    
    arma::vec g = rnorm();
    arma::vec w1 = arma::cumsum(r%g);
    if (!idx.is_empty()) {
      N = idx.n_elem;
      w1 = w1.elem(idx); 
    }
    arma::rowvec B = arma::sum(ic.each_col()%g, 0);
    arma::vec w2(N);
    if (!idx.is_empty()) {
      for (unsigned i=0; i<N; i++) {
	w2(i) = arma::as_scalar(B*eta.row(idx(i)).t());
      }
    } else {
      for (unsigned i=0; i<n; i++) {
	w2(i) = arma::as_scalar(B*eta.row(i).t());
      }
    }    
    return (w1+w2)/std::sqrt((double)n);
  }

  // Sample 'r' processes
  arma::mat cumres::sample(unsigned R, arma::uvec idx, bool quantiles) {
    qt.fill(0);    
    unsigned n = this->n;
    // arma::vec t0 = this->t;
    // if (!idx.is_empty()) {
    //   n = idx.n_elem;
    //   t0 = t0.elem(idx);
    // }
    arma::vec t0 = t;
    arma::mat qt(n, std::ceil(R*0.05));
    arma::mat res(R,2);
    for (unsigned i=0; i<R; i++) {
      arma::vec wi = this->sample(idx);
      res(i,0) = SupTest(wi);
      res(i,1) = L2Test(wi, t0);
      if (false) {
	for (unsigned j=0; j<n; j++) {
	  wi = abs(wi);
	  arma::rowvec qtj = qt.row(j);
	  unsigned i = qtj.index_min();
	  if (wi(j)>qt(j,i)) qt(j,i) = wi(j);
	}
      }
    }
    this->qt = qt;
    return res;
  }

}  // namespace target

