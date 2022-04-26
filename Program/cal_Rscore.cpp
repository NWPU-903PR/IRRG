#include <iostream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat cal_Rscore(GenericVector exp_ls,
               GenericVector back_ls,
               GenericVector adj_list,
               IntegerVector goi_vec) {
  int n_rows=exp_ls.size();
  int n_cols=5;
  double r_gy;
  double r_nogy;
  int out_num;
  int in_num;
  double total_mean;
  double part_mean;
  arma::mat R_info(n_rows,n_cols);
  for(int i =0;i<n_rows;i++){
    arma::sp_mat expgroup_temp = exp_ls[i];
    arma::mat expgroup = arma::mat(expgroup_temp);
    
    arma::sp_mat backgroup_temp = back_ls[i];
    arma::mat backgroup = arma::mat(backgroup_temp);
    
    arma::sp_mat adj_temp = adj_list[i];
    arma::mat adj = arma::mat(adj_temp);
    int goi_index  = goi_vec[i]-1;
    
    arma::rowvec in_deg = sum(adj,0);
    arma::colvec out_deg = sum(adj,1);
    
    out_num = out_deg[goi_index];
    in_num = in_deg[goi_index];
    if (backgroup.n_cols==0) {
      r_nogy=1;
      part_mean=arma::mean(expgroup.row(goi_index));
      total_mean=part_mean;
      R_info(i,0)=r_nogy;
      R_info(i,1)=out_num;
      R_info(i,2)=in_num;
      R_info(i,3)=total_mean;
      R_info(i,4)=part_mean;
      continue;

    }
    if (expgroup.n_cols==0) {
      r_gy=0;
      r_nogy=0;
      part_mean=0;
      total_mean=0;
      R_info(i,0)=r_gy;
      R_info(i,1)=out_num;
      R_info(i,2)=in_num;
      R_info(i,3)=total_mean;
      R_info(i,4)=part_mean;
      R_info(i,5)=r_nogy;
      continue;
    }
    
    arma::colvec bk_means=arma::mean(backgroup,1);
    arma::colvec exp_means=arma::mean(expgroup,1);
    arma::colvec fc_value = abs(bk_means-exp_means);
    arma::colvec r = fc_value;
    
    arma::rowvec in_deg_c=in_deg;
    in_deg_c=in_deg_c.replace(0,1);
    arma::mat M_c = adj;
    M_c.each_row()/= in_deg_c;
    
    arma::colvec f=fc_value/sum(fc_value);
    
    arma::rowvec up=arma::rowvec(size(in_deg));
    up.fill(1);
    arma::rowvec down=arma::rowvec(size(in_deg));
    down.fill(3);

    arma::colvec d = ((in_deg+up)/(in_deg+down)).t();
    arma::colvec one(in_deg.size(),arma::fill::ones);
    arma::colvec r_next = (one-d)%f+d%(M_c*r);
   
    double temp = sum(abs(r_next-r));
    while(temp>0.001){
      r_next = (one-d)%f+d%(M_c*r);
      temp = sum(abs(r_next-r));
      r=r_next;
      }
  
    r_nogy=r(goi_index);
    part_mean=arma::mean(expgroup.row(goi_index));
    total_mean=sum(expgroup.row(goi_index))/(expgroup.n_cols+backgroup.n_cols);
    R_info(i,0)=r_nogy;
    R_info(i,1)=out_num;
    R_info(i,2)=in_num;
    R_info(i,3)=total_mean;
    R_info(i,4)=part_mean;
    }
  return(R_info);
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//


