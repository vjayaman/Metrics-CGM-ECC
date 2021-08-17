#include <RcppArmadillo.h>
// [[ Rcpp :: depends ( RcppArmadillo )]]
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
using namespace Rcpp;



// // [[Rcpp::export]]
// vector<int> MultiplyVectorByScalar(vector<int> &v, int k){
//   transform(v.begin(), v.end(), v.begin(), [k](int &c){ return c*k; });
//   return v;
// }
// // from https://slaystudy.com/c-multiply-vector-by-scalar/
// // 1. We use std::transform to perform an operation on each element of the vector.
// // 2. The first 2 parameters, v.begin() and v.end() tells std::transform that we need to perform 
// // operations on each elements in range [ v.begin(), v.end() ).
// // 3. The third parameter is the initial iterator where we store/return the values after modification.
// // 4. The last parameter is the lambda function. The functions accept scalar k as a parameter. 
// // The function runs for each element of the vector. The function multiplies the vector element by k 
// // and returns the result.

/*** R
# v <- c(1, 2, 3, 4)
# k <- 4
# x <- MultiplyVectorByScalar(v, k)
# print(v)
# print(x)
*/

// [[Rcpp::export]]
double calculateEpi(int x, arma::mat dr_as, arma::mat epi_melt) {
  
  arma::colvec em_dr1 = epi_melt.col(0);
  arma::colvec em_dr2 = epi_melt.col(1);
  arma::colvec em_val = epi_melt.col(2);
  int dr;
  arma::uvec dr1_ids, dr2_ids; 
  
  arma::colvec clusters = dr_as.col(0); //cluster column of dr_assignments
  arma::colvec drvals = dr_as.col(1);
  arma::colvec nvals = dr_as.col(2);
  
  arma::uvec ids = find(clusters == x); // cluster x
  arma::colvec idrs = drvals.rows(ids);
  int nrows = idrs.n_rows;
  arma::colvec ins = nvals.rows(ids);
  
  for (int i = 0; i < nrows; i++) {
    dr = idrs[i];
    dr1_ids = find(em_dr1 == dr);
    dr2_ids = find(em_dr2 == dr);
    
    em_val.elem(dr1_ids) = em_val.elem(dr1_ids) * ins[i];
    em_val.elem(dr2_ids) = em_val.elem(dr2_ids) * ins[i];
  }
  
  double z = sum(em_val);
  return z;
}


// # cppFunction('List matchedDrs(DataFrame ccm, DataFrame dr_as) {
// #   int nrow = ccm.nrow();
// #   
// #   NumericVector x = ccm[0]; //cluster column of ccm
// #   NumericVector col1 = dr_as[0]; //cluster column of dr_assignments
// #   NumericVector col2 = dr_as[1]; //dr column of dr_assignments
// #   LogicalVector inds;
// #   
// #   for (int i = 0; i < nrow; i++) {
// #     inds = (col1 == x[i]); //rows in dr_assignments with cluster x[i]
// #   }
// #   return testv;
// # }')
// # 
// # inds <- matchedDrs(ccm, dr_assignments)
// 
// # NumericVector k = y[i]; //list of cluster members of cluster at row i
// 
// 
// 
// 
// # 
// # cppFunction('LogicalVector drsInCluster(DataFrame ccm, int i, DataFrame dr_as) {
// #   NumericVector x = ccm[0]; //cluster column of ccm
// #   
// #   NumericVector col1 = dr_as[0]; //cluster column of dr_assignments
// #   NumericVector col2 = dr_as[1]; //dr column of dr_assignments
// #   LogicalVector inds = (col1 == x[i]); //rows in dr_assignments with cluster x[i]
// #   
// #   return inds;
// # }')
// # 
// # cppFunction('double sumEpiVals(DataFrame epi_melt) {
// #   NumericVector value = epi_melt[2];
// #   NumericVector n1 = epi_melt[3];
// #   NumericVector n2 = epi_melt[4];
// #   
// #   NumericVector value2 = value * n1 * n2;
// #   return sum(value2);
// # }')
// # 
// # cppFunction('double calculateEpi(int x, arma::mat dr_as, arma::mat epi_melt) {
// # 
// #   arma::colvec em_dr1 = epi_melt.col(0);
// #   arma::colvec em_dr2 = epi_melt.col(1);
// #   arma::colvec em_val = epi_melt.col(2);
// #   int dr;
// #   arma::uvec dr1_ids, dr2_ids; 
// #   
// #   arma::colvec clusters = dr_as.col(0); //cluster column of dr_assignments
// #   arma::colvec drvals = dr_as.col(1);
// #   arma::colvec nvals = dr_as.col(2);
// #   
// #   arma::uvec ids = find(clusters == x); // cluster x
// #   arma::colvec idrs = drvals.rows(ids);
// #   int nrows = idrs.n_rows;
// #   arma::colvec ins = nvals.rows(ids);
// #   
// #   for (int i = 0; i < nrows; i++) {
// #     dr = idrs[i];
// #     dr1_ids = find(em_dr1 == dr);
// #     dr2_ids = find(em_dr2 == dr);
// #     
// #     em_val.elem(dr1_ids) = em_val.elem(dr1_ids) * ins[i];
// #     em_val.elem(dr2_ids) = em_val.elem(dr2_ids) * ins[i];
// #   }
// #   
// #   double z = sum(em_val);
// #   return z;
// # }', depends = "RcppArmadillo")
// 
// # x1 <- calculateEpi(ccm, 2, dr_assignments, epi_melt)
// # x2 <- calculateEpi(ccm, 2, dr_assignments, epi_melt)
