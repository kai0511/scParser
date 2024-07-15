// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/scParser_types.h"
#include <iostream>
#include "utils.h"

using std::pow;

field<uvec> generate_batches(const unsigned int& sample_size, const unsigned int& num_batch){
    field<uvec> batches(num_batch);
    uvec shuffled = randperm(sample_size);
    unsigned int batch_size = sample_size/num_batch;
    for(unsigned int i = 0; i < num_batch; i++){
        if(i == num_batch - 1){
            batches(i) = shuffled(span(i * batch_size, sample_size - 1));
        }else {
            batches(i) = shuffled(span(i * batch_size, (i + 1) * batch_size - 1));
        }
    }
    return batches;
}

double objective(const mat& X, const vec& y, const vec& beta,
                 const double& lambda, const double& alpha){
    vec residual = y - X * beta;
    double loss = sum(square(residual))/2;
    double l2l1_reg = (1 - alpha) * lambda * sum(square(beta))/2 + alpha * lambda * sum(abs(beta));
    loss += l2l1_reg;
    return loss;
}

double compute_loss(const vec& residual, const vec& beta, const double& lambda, const double& alpha){
    double loss = sum(square(residual))/2 + (1 - alpha) * lambda * sum(square(beta))/2 + alpha * lambda * sum(abs(beta));
    return loss;
}

void predict(const mat& row_factor, const mat& column_factor, const mat& cell_factor, const mat& gene_factor, mat& predictions) {
    predictions = row_factor * column_factor + cell_factor * gene_factor;
}

void evaluate(mat& residual, const uvec& train_idx, const uvec& test_idx, 
              double& sum_residual, double& train_rmse, double& test_rmse, const int& tuning, 
              const int& iter = 0, const int& verbose = 1){

    // residual = data - predictions;
    if(tuning == 0){
        sum_residual = accu(sum(square(residual), 1));
        train_rmse = std::sqrt(sum_residual/(residual.n_cols * residual.n_rows));
    }else{
        sum_residual = sum(square(residual.elem(train_idx)));
        train_rmse = std::sqrt(sum_residual/train_idx.n_elem);
        test_rmse = std::sqrt(mean(square(residual.elem(test_idx))));
    } 

    if (verbose == 1){
        cout << "SC2 iter " << iter << ": train rmse = " << train_rmse << endl;

        if(tuning == 1){
            cout << "SC2 iter " << iter << ": test rmse = " << test_rmse << endl;
        }
    }
}

double compute_loss(const field<mat>& cfd_factor, const mat& column_factor, const mat& cell_factor, const mat& gene_factor, 
                    const double& lambda1, const double& lambda2, const double& alpha, double& sum_residual, const int& verbose = 1){

    // l2 penalty 
    double row_reg = 0.0;
    for(unsigned int i = 0; i < cfd_factor.n_elem; i++){
        row_reg += lambda1 * pow(norm(cfd_factor(i), "F"), 2);
    }
    double col_reg = lambda1 * pow(norm(column_factor, "F"), 2);
    
    double cell_l1_reg = lambda2 * alpha * sum(sum(abs(cell_factor), 1));   //  l1 penalty
    double cell_l2_reg = lambda2 * (1 - alpha) * pow(norm(cell_factor, "F"), 2); //  l2 penalty

    double loss = sum_residual/2 + row_reg/2 + col_reg/2 + cell_l1_reg + cell_l2_reg/2; 

    if(verbose == 1) {
        cout << "total_residual" << '\t' << sum_residual/2 << ";" << '\n'
             << "row_reg_loss:" << '\t' << row_reg/2 << ";" << '\n'
             << "col_reg_loss:" << '\t' << col_reg/2 << ";" << '\n'
             << "cell_l2_reg:" << '\t' << cell_l2_reg/2 << ";" << '\n'
             << "l1_reg_loss:" << '\t' << cell_l1_reg << "." << endl;
    }
    return loss;
}
