#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
int clusterSizes(DataFrame x, int h) {
  NumericVector cluster = x[0];
  NumericVector y = x[2];
  int nrow = x.nrow();
  int vals = 0;
  
  for (int i = 0; i < nrow; i++) {
    if (cluster[i] == h) {
      vals = vals + y[i];
    }
  }
  return vals;
}


/*** R
# libs <- c("R6","testit","optparse","magrittr","dplyr","tibble","readr",
#           "reshape2","fossil","tidyr","purrr", "data.table", "Rcpp", "RcppArmadillo")
# y <- lapply(libs, require, character.only = TRUE)
# assert("All packages loaded correctly", all(unlist(y)))
# sourceCpp("scripts/epicohversions.cpp")
# source("scripts/ECC/ecc_functions.R")

# g_cuts <- readRDS("intermediate_data/g_cuts.Rds")
# epi_melt <- readRDS("intermediate_data/epi_melt.Rds")
# t1 <- Sys.time()
# df1 <- epiCohesionV1(g_cuts, epi_melt)
# t2 <- Sys.time()
# df2 <- epiCohesionV2(g_cuts, epi_melt)
# t3 <- Sys.time()
# 
# paste0("epiCohV1 took ", t2 - t1)
# paste0("epiCohV2 took ", t3 - t2)
*/

// [[Rcpp::export]]
NumericVector calculateEpiV2(NumericVector clusterlist, arma::mat dr_as, arma::mat epi_melt) {
  
  int cl_len = clusterlist.size();
  NumericVector y(cl_len);
  
  for (int j = 0; j < cl_len; j++) {
    int x = clusterlist[j];
    
    arma::colvec em_dr1 = epi_melt.col(0);
    arma::colvec em_dr2 = epi_melt.col(1);
    arma::colvec em_val = epi_melt.col(2);
    
    arma::colvec em_flag1 = epi_melt.col(3);
    arma::colvec em_flag2 = epi_melt.col(4);
    int dr;
    arma::uvec dr1_ids, dr2_ids; 
    
    arma::colvec clusters = dr_as.col(0); //cluster column of dr_assignments
    arma::colvec drvals = dr_as.col(1);
    arma::colvec nvals = dr_as.col(2);
    
    arma::uvec ids = find(clusters == x); // cluster x
    arma::colvec idrs = drvals.elem(ids);
    int nrows = idrs.n_rows;
    arma::colvec ins = nvals.elem(ids);
    
    for (int i = 0; i < nrows; i++) {
      dr = idrs[i];
      dr1_ids = find(em_dr1 == dr);
      dr2_ids = find(em_dr2 == dr);
      
      em_val.elem(dr1_ids) = em_val.elem(dr1_ids) * ins[i];
      em_val.elem(dr2_ids) = em_val.elem(dr2_ids) * ins[i];
      em_flag1.elem(dr1_ids).fill(1);
      em_flag2.elem(dr2_ids).fill(1);
    }
    
    arma::colvec flags = em_flag1 + em_flag2;
    arma::uvec final_ids = find(flags == 2);
    arma::colvec actual = em_val.elem(final_ids);
    
    em_flag1.fill(0);
    em_flag2.fill(0);
    y[j] = sum(actual);
  }
  
  return y;
}

// [[Rcpp::export]]
double allClusters(arma::colvec clusters) {
  
  for (int x : clusters) {
    Rprintf("%d \n", x);
  }
  return 0;
}

// [[Rcpp::export]]
NumericMatrix transformTempDists(NumericMatrix dm, double min_dist, double max_dist) {
  double lmin, lmax;
  int nrows = dm.nrow();
  lmin = log10(min_dist + 10);
  lmax = log10(max_dist + 10);

  dm = dm + 10;

  for (int j = 0; j < nrows; j++) {dm.row(j) = log10(dm.row(j));}

  if (lmax == 0) {
    dm = 0;
  }else {
    for (int j = 0; j < nrows; j++) {
      dm.row(j) = (dm.row(j) - lmin) / (lmax - lmin);
    }
  }
  return dm;
}


// [[Rcpp::export]]
NumericMatrix transformGeoDists(NumericMatrix dm, double min_dist, double max_dist) {
  double lmin, lmax;
  int nrows = dm.nrow();
  lmin = log10(min_dist + 10);
  lmax = log10(max_dist + 10);

  dm = dm + 10;

  for (int j = 0; j < nrows; j++) {dm.row(j) = log10(dm.row(j));}

  if (lmax == 1) {
    dm = 0;
  }else {
    for (int j = 0; j < nrows; j++) {
      dm.row(j) = (dm.row(j) - lmin) / (lmax - lmin);
    }
  }
  return dm;
}

// [[Rcpp::export]]
DataFrame epiMat(NumericMatrix temp, NumericMatrix geo, double tau, double gamma) {
  int dim = temp.nrow(); // should be same for all dims, by how they were made
  NumericMatrix em(dim, dim);
  
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      em(i,j) = 1 - sqrt(pow(temp(i,j), 2)*tau + pow(geo(i,j), 2)*gamma);
    }
  }
  return em;
}