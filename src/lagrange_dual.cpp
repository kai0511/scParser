// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/scParser_types.h"
#include <iostream>
#include <omp.h>
#include "lagrange_dual.h"

using namespace arma;

// [[Rcpp::export]]
mat lagrange_dual(const mat& gram, const mat& cdx, const int& reg = 1, 
                  const int& max_iter = 100, const double& tol = 1e-4){
    
    vec lambda = ones(gram.n_cols);
    vec gradient, new_lambda;
    mat lmbd, inv_gram, base, hassian;
    
    // mat gram = trans(codes) * codes;
    // mat cdx = trans(codes) * X;
    for(int iter = 0; iter < max_iter; iter++){
        lmbd = diagmat(lambda);
        inv_gram = inv(gram + lmbd);
        base = solve(gram + lmbd, cdx, solve_opts::likely_sympd);
        
        // pre_dual_loss = dual_loss;
        // dual_loss = objective(X, codes, base, lmbd);
        // cout << "Delta dual loss:" << dual_loss - pre_dual_loss << endl;
        
        // gradient
        gradient = (sum(square(base), 1) - reg)/2;

        // hessian
        hassian = - base * trans(base) % inv_gram;

        new_lambda = lambda - solve(hassian, gradient, solve_opts::likely_sympd);

        if(sum(square(lambda - new_lambda)) < tol){
            break;
        }
        lambda = new_lambda;
    }
    base = solve(gram + diagmat(new_lambda), cdx, solve_opts::likely_sympd);
    return base;
}
