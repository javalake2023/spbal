// spbal.cpp

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]

#include <Rcpp.h>
#include <RcppThread.h>

using Rcpp::IntegerVector;
using Rcpp::NumericVector;
using namespace Rcpp;

#include <cmath>

// not used yet, but here when we need them.
#define RCPPTHREAD_OVERRIDE_COUT 1   // std::cout override
#define RCPPTHREAD_OVERRIDE_THREAD 1 // std::thread override


//' @name log_a_to_base_b
//'
//' @title Compute the log of a to base b.
//'
//' @description Compute the log of a to base b.
//'
//' @details This function was written by Phil Davies.
//'
//' @param a Integer to find the log to base b of.
//' @param b Base
//'
//' @return The log of a to base b.
//'
//' @keywords internal
double log_a_to_base_b(long long a, int b){

  return std::log2(a) / std::log2(b);
}


//' @name mod
//'
//' @title Need to fill in.
//'
//' @description Need to fill in.
//'
//' @details This function was written by Phil Davies.
//'
//' @param a Need to fill in.
//' @param n Need to fill in.
//'
//' @return Need to fill in.
//'
//' @keywords internal
template<typename T>
T mod(T a, int n){
  return a - floor(a / n) * n;
}


//' @name removeDuplicates
//'
//' @title Need to fill in.
//'
//' @description Need to fill in.
//'
//' @details This function was written by Phil Davies.
//'
//' @param vec Need to fill in.
//'
//' @return Need to fill in.
//'
//' @keywords internal
Rcpp::NumericVector removeDuplicates(Rcpp::NumericVector vec){
  // sort vector
  std::sort(vec.begin(), vec.end());
  // remove duplicates
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
  return vec;
}


Rcpp::IntegerVector sample_int_v2(int max_int, int num_seeds) {
  Rcpp::IntegerVector v = Rcpp::sample(max_int, num_seeds);
  return v;
}


//' @name cppBASpts
//'
//' @title Generate numbers from a Halton Sequence.
//'
//' @description For efficiency, this function can generate points along a random start
//' Halton Sequence for a predefined Halton.
//'
//' @details This function was first written in R by Blair Robertson, subsequently it was
//' re-written in C/C++ by Phil Davies.
//'
//' @param n Number of points required.
//' @param seeds Random starting point in each dimension.
//' @param bases Co-prime base for the Halton Sequence.
//' @param verbose A boolean indicating whether informational messages are to be issued.
//'
//' @return Matrix with the columns, order of points, x in [0,1) and y in [0,1)
//'
//' @examples
//' # First 10 points in the Halton Sequence for base 2,3
//' spbal::cppBASpts(n = 10)
//' # First 10 points in the Halton Sequence for base 2,3 with
//' # starting point at the 15th and 22nd index.
//' spbal::cppBASpts(n = 10, seeds = c(14, 21))
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cppBASpts(int n = 10,
                     IntegerVector seeds = Rcpp::IntegerVector::create(),
                     NumericVector bases = Rcpp::NumericVector::create(),
                     bool verbose = false)
{
  // set default seeds values
  if (seeds.size() == 0){
    //seeds = sample_int(2, 0, 62208);
    seeds = sample_int_v2(62208, 2);
  }
  // set default bases values
  if (bases.size() == 0){
    bases = {2, 3};
  }

  //if(verbose){
  //  RcppThread::Rcout << "cppBASpts() n     : " << n     << std::endl;
  //  RcppThread::Rcout << "cppBASpts() seeds : " << seeds << std::endl;
  //  RcppThread::Rcout << "cppBASpts() bases : " << bases << std::endl;
  //}

  // initialise variables
  int d = bases.length();
  int u, b;

  Rcpp::NumericVector xk;
  Rcpp::List          xklist;
  Rcpp::NumericMatrix pts(n, d);

  if (seeds.length() != d){
    //if(verbose){
    //  RcppThread::Rcout << "cppBASpts() seeds.length() != d : " << std::endl;
    //}
    seeds = Rcpp::rep(seeds[1], d);
  }

  //########### Main Loop #########################################
  for (int i = 0; i < d; i++){
    b = bases[i];
    u = seeds[i];
    // k <- u:(u+n-1);
    Rcpp::IntegerVector ik = Rcpp::seq(u, (u+n-1));
    Rcpp::NumericVector k = as<Rcpp::NumericVector>(ik);
    xk = mod(k, b) / b;
    xklist.push_back(removeDuplicates(clone(xk)));

    for (int j = 0; j < (std::ceil(log_a_to_base_b(u + n, b)) + 2); j++){
      Rcpp::NumericVector tmp1 = Rcpp::floor(k / std::pow(b, j + 1));
      //int tmp000 = std::pow(b, (j + 2));
      NumericVector tmp2 = mod(tmp1, b) / std::pow(b, (j + 2));
      xk = xk + tmp2;
    } // end for j

    // point to column i
    NumericMatrix::Column thisCol = pts(_, i);
    // propagate changes to pts.
    thisCol = thisCol + xk;

  } // end for i

  //return pts;
  return Rcpp::List::create(_["pts"]    = pts,
                            _["xklist"] = xklist,
                            _["seeds"]  = seeds);
}


//' @name cppRSHalton_br
//'
//' @title Generate numbers from a Halton Sequence with a random start
//'
//' @description For efficiency, this function can generate points along a random start
//' Halton Sequence for a predefined Halton.
//'
//' @details This function was first written in R by Paul van Dam-Bates for the
//' package BASMasterSample. Subsequently it was written in C/C++ by Phil Davies.
//'
//' @param n Number of points required
//' @param bases Co-prime base for the Halton Sequence
//' @param seeds Random starting point in each dimension
//' @param verbose A boolean indicating whether informational messages are to be issued.
//'
//' @return Matrix with the columns, order of point, x in [0,1) and y in [0,1).
//'
//' @examples
//' # First 10 points in the Halton Sequence for base 2,3
//'  spbal::cppRSHalton_br(n = 10)
//' # First 10 points in the Halton Sequence for base 2,3 with
//' # starting point at the 15th and 22nd index.
//'  spbal::cppRSHalton_br(n = 10, seeds = c(14, 21))
//'
//' @export
// [[Rcpp::export(rng = false)]]
Rcpp::List cppRSHalton_br(int n = 10,
                          NumericVector bases = Rcpp::NumericVector::create(),
                          NumericVector seeds = Rcpp::NumericVector::create(),
                          bool verbose = false)
{
  // defaults: n = 10, bases = c(2, 3), seeds = c(0, 0)

  //if (seeds.size() == 0){
  //  seeds = {0, 0};
  //}
  if (bases.size() == 0){
    bases = {2, 3};
  }

  //if(verbose){
  //  RcppThread::Rcout << "cppRSHalton_br() n     : " << n     << std::endl;
  //  RcppThread::Rcout << "cppRSHalton_br() seeds : " << seeds << std::endl;
  //  RcppThread::Rcout << "cppRSHalton_br() bases : " << bases << std::endl;
  //}

  // initialize variables.
  long long UpLim = std::pow(10, 15);
  long long m = 0;
  long long tmpm, u2;

  int d = bases.length();
  int b;

  double min = 0.0; // Minimum value for runif
  double max = 1.0; // Maximum value for runif

  std::vector<double> u;

  Rcpp::NumericVector xk;
  Rcpp::NumericMatrix pts(n, d);
  Rcpp::NumericVector k(n);
  Rcpp::List          xklist;

  if (seeds.size() == 0){
    // seeds = numeric(d)
    seeds[d];
    // u = runif(d)
    // Generate random values using uniform distribution
    for (int i = 0; i < d; i++) {
      // Generate random uniform values (respects the current R set.seed())
      //RcppThread::Rcout << "cppRSHalton_br() u : " << R::runif(min, max) << std::endl;
      u.push_back(R::runif(min, max));
    } // end for i

    // for each element of bases
    for (int i = 0; i < d; i++) {
      b = bases[i];
      m = 0;
      int j = 1;
      // while (m + (b-1)*(b^(j-1)) <= UpLim) {
      while (m + (b-1) * std::pow(b, j-1) <= UpLim){
        // m = m + (floor(u[i]*(b^j)) %% b)*(b^(j-1))
        tmpm = std::floor(u[i] * std::pow(b, j));
        m = m + (mod(tmpm, b) * std::pow(b, j-1));
        j++;
      }
      seeds[i] = m;
    } // end for i
  } // end if seeds.size()

  //########### Main Loop #########################################
  for (int i = 0; i < d; i++){
    b = bases[i];
    u2 = seeds[i];

    // Populate the vector with values from u to (u+n-1), R: k <- u:(u+n-1);
    for (int i = 0; i < n; i++) {
      k[i] = (u2 + i);
    }

    xk = mod(k, b) / b;
    xklist.push_back(removeDuplicates(clone(xk)));

    for (int j = 1; j <= (std::ceil(log_a_to_base_b(u2 + n, b)) + 2); j++){
      //
      NumericVector tmp1 = Rcpp::floor(k / std::pow(b, j));
      NumericVector tmp2 = mod(tmp1, b) / std::pow(b, (j + 1));
      xk = xk + tmp2;
    } // end for j

    // point to column i
    NumericMatrix::Column thisCol = pts(_, i);
    // propagate changes to pts.
    thisCol = thisCol + xk;
  } // end for i

  //return pts;
  return Rcpp::List::create(_["pts"]    = pts,
                            _["xklist"] = xklist,
                            _["seeds"]  = seeds);
}
