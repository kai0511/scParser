// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/SR2_types.h"
#include <iostream>

vec nnls(const mat &X, const subview_col<double> &y, const double lambda, int max_iter = 500, double tol = 1e-6, bool verbose = false){
    /*
     * Description: sequential Coordinate-wise algorithm for non-negative least square regression A x = b, s.t. x >= 0
     *  Reference: http://cmp.felk.cvut.cz/ftp/articles/franc/Franc-TR-2005-06.pdf 
     */

    vec mu = -X.t() * y;
    mat H = X.t() * X;
    vec beta(X.n_cols), beta0(X.n_cols);
    beta.fill(0);
    beta0.fill(-9999);

    int i = 0;
    double tmp, beta_sum = 0.0;

    while(i < max_iter && max(abs(beta - beta0)) > tol) {
        beta0 = beta;

        for (int k = 0; k < X.n_cols; k++) {

            // tmp = beta[k] - mu[k] / H.at(k,k);
            tmp = beta[k] - (mu[k] + lambda * beta_sum)/(H.at(k,k) + lambda)

            if (tmp < 0) tmp = 0;

            if (tmp != beta[k]) {
                mu += (tmp - beta[k]) * H.col(k);
                beta_sum += tmp - beta[k];
            }
            beta[k] = tmp;
        }
        ++i;
    }
    return beta;
}