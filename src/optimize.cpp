// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins("cpp11")]]

#include "../inst/include/SR2_types.h"
#include <iostream>
#include <omp.h>
#include "lagrange_dual.h"
#include "coordinate_descent.h"
#include "utils.h"
#include <ctime>

using Rcpp::List;
using Rcpp::Named;

void tune_cfd_row(const mat& residual, const umat& indicator, mat& updating_factor, const mat& c_factor, 
                  const uvec& updating_confd, const mat& gram, const double& lambda, const unsigned int& n_cores = 10) {
    /*
        fix column parameters and update row factors 
    args:
        @dimension must be in (1, 2), denoting which row_factor will be updated
        @n_cores number of cores for parallel computation
    */
    
    uvec levels = unique(updating_confd);

    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
    #endif
    for(unsigned int i = 0; i < levels.size(); i++) {

        uvec non_zeros, zero_idx;
        mat XtX = zeros(c_factor.n_rows, c_factor.n_rows);
        vec outcome, Xty = zeros<vec>(c_factor.n_rows);

        uvec ids = find(updating_confd == levels(i));
        for(unsigned int k = 0; k < ids.n_elem; k++){
            
            non_zeros = find(trans(indicator.row(ids(k))));   
            zero_idx = find(trans(indicator.row(ids(k))) == 0);           
            outcome = trans(residual.row(ids(k)));
            outcome = outcome(non_zeros);

            // XtX += c_factor.cols(non_zeros) * trans(c_factor.cols(non_zeros));
            XtX += gram - c_factor.cols(zero_idx) * trans(c_factor.cols(zero_idx));
            Xty += c_factor.cols(non_zeros) * outcome;
            // Xty += Xtys.col(ids(k)) - c_factor.cols(zero_idx) * outcome;
        }
        
        XtX.diag() += lambda;
        updating_factor.row(i) = trans(solve(XtX, Xty, solve_opts::likely_sympd));
    }
}

void fit_cfd_row(const mat& residual, mat& updating_factor, const mat& c_factor, 
                 const mat& updating_confd, const vec& cfd_cnts, const mat& gram, const double& lambda, const unsigned int& n_cores = 5) {
    /*
        fix column parameters and update row factors 
    args:
        @dimension must be in (1, 2), denoting which row_factor will be updated
        @n_cores number of cores for parallel computation
    */

    unsigned int nlevels = updating_confd.n_cols;
    mat Xtys = c_factor * (residual.t() * updating_confd);
    
    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 1)
    #endif
    for(unsigned int i = 0; i < nlevels; i++) {
        // vec Xty = Xtys * updating_confd.col(i);
        updating_factor.row(i) = trans(solve(cfd_cnts(i) * gram + lambda * eye(size(gram)), Xtys.col(i), solve_opts::likely_sympd));
    }
}

void optimize_cfd_col(const mat& residual, const umat& indicator, const mat& row_factor, mat& c_factor, 
                      const double& lambda, const int tuning, const int n_cores = 10){
    
    if(tuning == 1){
        mat gram = trans(row_factor) * row_factor;

        #if defined(_OPENMP)
            #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 10)
        #endif
        for(unsigned int i = 0; i < residual.n_cols; i++) {
            mat XtX;
            uvec row_selected = find(indicator.col(i));

            if(2*row_selected.n_elem  <= indicator.n_rows){
                // XtX = sum(feature_space.slices(row_selected), 2);
                XtX = trans(row_factor.rows(row_selected)) * row_factor.rows(row_selected);
            }else{
                uvec unselected = find(indicator.col(i) == 0);
                // XtX = sum(feature_space.slices(unselected), 2);
                XtX = gram - row_factor.rows(unselected).t() * row_factor.rows(unselected);
            }
            XtX.diag() += lambda; 

            vec outcome = residual.col(i);
            outcome = outcome(row_selected);
            c_factor.col(i) = solve(XtX, trans(row_factor.rows(row_selected)) * outcome, solve_opts::likely_sympd);
        }

    } else if(tuning == 0) {

        mat XtX = row_factor.t() * row_factor;
        mat Xty = row_factor.t() * residual;
        XtX.diag() += lambda; 
        mat inv_XtX = inv(XtX);

        #if defined(_OPENMP)
            #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 50)
        #endif
        for(unsigned int i = 0; i < residual.n_cols; i++) {
            // c_factor.col(i) = solve(XtX, Xty.col(i), solve_opts::likely_sympd);
            c_factor.col(i) = inv_XtX * Xty.col(i);
        }
    } else {
        cout << "Parameter tuning should be either 0 or 1!" << endl;
        exit(1);
    }
}

// [[Rcpp::export]]
void tune_coding(const mat& residual, const umat& indicator, mat& row_factor, const mat& c_factor, 
                     const double lambda, const double alpha, int tuning = 1, int n_cores = 20, double tol = 1e-5){
    
    mat tc_factor = trans(c_factor);

    mat gram = c_factor * tc_factor;
    // cube feature_space(c_factor.n_rows, c_factor.n_rows, c_factor.n_cols);
    // for(unsigned int i = 0; i < c_factor.n_cols; i++){
    //     feature_space.slice(i) = c_factor.col(i) * tc_factor.row(i);
    // }

    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 100)
    #endif
    for(unsigned int i = 0; i < residual.n_rows; i++) {
        mat XtX;
        uvec row_selected = find(trans(indicator.row(i)));

        if(row_selected.n_elem*2 <= indicator.n_cols){
            // XtX = sum(feature_space.slices(row_selected), 2);
            XtX = c_factor.cols(row_selected) * tc_factor.rows(row_selected);
        }else{
            uvec unselected = find(trans(indicator.row(i)) == 0);
            // XtX = sum(feature_space.slices(unselected), 2);
            XtX = gram - c_factor.cols(unselected) * tc_factor.rows(unselected);
        }
        vec outcome = trans(residual.row(i));
        outcome = outcome(row_selected);
        // vec Xty = c_factor.cols(row_selected) * outcome;

        if(alpha == 0){
            XtX.diag() += lambda;
            row_factor.row(i) = trans(solve(XtX, c_factor.cols(row_selected) * outcome, solve_opts::likely_sympd));
        } else {
            row_factor.row(i) = trans(strong_coordinate_descent(tc_factor.rows(row_selected), outcome, trans(row_factor.row(i)), lambda, alpha, XtX, c_factor.cols(row_selected) * outcome, tol));
        }
    }
}

// [[Rcpp::export]]
void fit_coding(const mat& residual, mat& row_factor, const mat& c_factor, 
                const double lambda, const double alpha, int tuning = 1, int n_cores = 20, double tol = 1e-5){

    mat tc_factor = trans(c_factor); 
    mat XtX = c_factor * tc_factor;

    #if defined(_OPENMP)
        #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 100)
    #endif
    for(unsigned int i = 0; i < residual.n_rows; i++) {
        // vec Xty = c_factor * trans(residual.row(i));

        if(alpha == 0){
            XtX.diag() += lambda;
            row_factor.row(i) = trans(solve(XtX, c_factor * trans(residual.row(i)), solve_opts::likely_sympd));
        } else {
            vec outcome = trans(residual.row(i));
            row_factor.row(i) = trans(strong_coordinate_descent(tc_factor, outcome, trans(row_factor.row(i)), lambda, alpha, XtX, c_factor * trans(residual.row(i)), tol));
        }
    }
}

// [[Rcpp::export]]
void optimize_base(const mat& data, const umat& indicator, const mat& row_factor, mat& c_factor, 
                   const unsigned int& tuning, const double& tol = 1e-5, const unsigned int& iter = 10){
    
    if(tuning == 0){
        mat gram = trans(row_factor) * row_factor;
        mat cdx = trans(row_factor) * data;
        // c_factor = lagrange_dual(data, row_factor, 1, 100, 1e-4);
        c_factor = lagrange_dual(gram, cdx, 1, 100, 1e-4);

    }else{
        // if tuning, then block coordinate descent is used to optimize the problem.
        mat residual = data - row_factor * c_factor;
        uvec zero_idx = find(indicator == 0);
        residual.elem(zero_idx).zeros();

        double pre_loss, loss = accu(pow(residual, 2));
        // rowvec feature, XtX, Xty;

        // mat tran_row_factor = trans(row_factor);
        mat XtX = trans(pow(row_factor, 2)) * indicator;
        
        for(unsigned int i = 0; i < iter; i++) {
            
            for(unsigned int k = 0; k < c_factor.n_rows; k++) {
                residual += row_factor.col(k) * c_factor.row(k);
                residual.elem(zero_idx).zeros();
                
                c_factor.row(k) = normalise(trans(row_factor.col(k)) * residual / XtX.row(k));
                residual -= row_factor.col(k) * c_factor.row(k); 
            }

            pre_loss = loss;
            residual.elem(zero_idx).zeros();
            loss = accu(pow(residual, 2));

            if(abs(pre_loss - loss) < tol){
                break;
            }
        }
    }
}

// [[Rcpp::export]]
List partial_optimize(const mat& data, const umat& train_indicator, List& cfd_factors, mat column_factor, const umat& cfd_indicators, 
              const double lambda1 = 0.1, const int tuning = 1, const double global_tol = 1e-10, const unsigned int max_iter = 10000){
    
    cout.precision(12);
    // time_t tic, toc; 
    unsigned int i, iter = 0, cfd_num = cfd_factors.size();
    uvec train_idx, test_idx;
    double loss = 0.0, pre_loss, delta_loss, sum_residual, train_rmse, test_rmse; 
    mat gram, residual, sub_pred, row_factor = zeros(data.n_rows, column_factor.n_rows);

    // check whether the number of the confounding matrices is equal to the number of confounding indicators.
    if(cfd_num != cfd_indicators.n_cols){
        cout << "The number of confounding matrices should be the same the number of indicators." << endl;
        exit(1);
    }

    // place indices of confounders into arma::field for computational consideration
    field<mat> index_matrices(cfd_num);
    field<vec> confd_counts(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        uvec levels = unique(cfd_indicators.col(i));
        mat cfd_idx = zeros(cfd_indicators.n_rows, levels.n_elem);

        for(unsigned int k = 0; k < levels.n_elem; k++) {
            vec tmp = cfd_idx.col(k);
            tmp.elem(find(cfd_indicators.col(i) == levels(k))).ones();
            cfd_idx.col(k) = tmp;
        }
        index_matrices(i) = cfd_idx;
        confd_counts(i) = trans(sum(cfd_idx));
    }

    // move confounding matrices from Rcpp::List into arma::field
    field<mat> cfd_matrices(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        Rcpp::NumericMatrix temp = cfd_factors[i];
        cfd_matrices(i) = mat(temp.begin(), temp.nrow(), temp.ncol(), false);
        row_factor += index_matrices(i) * cfd_matrices(i);
    }

    // find the indices for training and testing elements 
    train_idx = find(train_indicator);
    test_idx = find(train_indicator == 0);

    // check the fitting with initial values
    residual = data - row_factor * column_factor;
    evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);

    // compute loss
    loss = sum_residual/2;
    for(unsigned int i = 0; i < cfd_matrices.n_elem; i++){
        loss += lambda1 * pow(norm(cfd_matrices(i), "F"), 2)/2;
    }
    loss += lambda1 * pow(norm(column_factor, "F"), 2)/2;

    while(iter <= max_iter) {

        if(iter % 10 == 0){
            cout << "Iteration " << iter << " ---------------------------------" << endl;
        }
        
        // update all confonding matrices
        gram = column_factor * column_factor.t();
        for(i = 0; i < cfd_num; i++){

            residual += index_matrices(i) * cfd_matrices(i) * column_factor;
            if(tuning == 1){
                tune_cfd_row(residual, train_indicator, cfd_matrices(i), column_factor, cfd_indicators.col(i), gram, lambda1, tuning);
            }else if(tuning == 0){
                fit_cfd_row(residual, cfd_matrices(i), column_factor, index_matrices(i), confd_counts(i), gram, lambda1, tuning);
            }else{
                cout << "Parameter tuning should be either 0 or 1!" << endl;
                exit(1);
            }

            // a trick to reduce computational burden
            if(i != cfd_num - 1){
                residual -= index_matrices(i) * cfd_matrices(i) * column_factor;
            }

	        evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
        }

        // update columm_factor
        row_factor.zeros();
        for(i = 0; i < cfd_num; i ++) {
            row_factor += index_matrices(i) * cfd_matrices(i);
        }

        optimize_cfd_col(data, train_indicator, row_factor, column_factor, lambda1, tuning, 30);
        residual = data - row_factor * column_factor;

        // check fitting every 10 steps
        if(iter % 10 == 0){
            evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);

            // compute loss
            pre_loss = loss;
            loss = sum_residual/2;
            for(unsigned int i = 0; i < cfd_matrices.n_elem; i++){
                loss += lambda1 * pow(norm(cfd_matrices(i), "F"), 2)/2;
            }
            loss += lambda1 * pow(norm(column_factor, "F"), 2)/2; 

            delta_loss = pre_loss - loss;
            cout << "Delta loss for iter " << iter << ":" << delta_loss << endl;

            if(delta_loss/pre_loss < global_tol){
                break;
            }
        }
        iter++;
    }

    // put the updated confounding matrices into a R list
    List row_matrices;
    for(i = 0; i < cfd_num; i ++) {
        row_matrices["factor" + std::to_string(i)] = cfd_matrices(i);
    }

    return List::create(Named("train_indicator") = train_indicator,
                        Named("row_matrices") = row_matrices,
                        Named("column_factor") = column_factor, 
                        Named("test_rmse") = test_rmse, 
                        Named("loss") = loss);
}

// [[Rcpp::export]]
List optimize(const mat& data, const umat& train_indicator, List& cfd_factors, mat column_factor, const umat& cfd_indicators, 
              mat cell_factor, mat gene_factor, const double lambda1 = 0.1, const double lambda2 = 0.01,
              const double alpha = 1.0, const int tuning = 1, const double global_tol = 1e-10, const double sub_tol = 1e-5, const unsigned int max_iter = 10000){
    
    cout.precision(12);
    // time_t tic, toc; 
    unsigned int i, iter = 0, cfd_num = cfd_factors.size();
    uvec train_idx, test_idx;
    double loss, pre_loss, delta_loss, sum_residual, train_rmse, test_rmse; // decay = 1.0; 
    mat gram, residual, sub_pred, row_factor = zeros(data.n_rows, column_factor.n_rows) , predictions = zeros(size(data));

    // check whether the number of the confounding matrices is equal to the number of confounding indicators.
    if(cfd_num != cfd_indicators.n_cols){
        cout << "The number of confounding matrices should be the same the number of indicators." << endl;
        exit(1);
    }

    // place indices of confounders into arma::field for computational consideration
    field<mat> index_matrices(cfd_num);
    field<vec> confd_counts(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        uvec levels = unique(cfd_indicators.col(i));
        mat cfd_idx = zeros(cfd_indicators.n_rows, levels.n_elem);

        for(unsigned int k = 0; k < levels.n_elem; k++) {
            vec tmp = cfd_idx.col(k);
            tmp.elem(find(cfd_indicators.col(i) == levels(k))).ones();
            cfd_idx.col(k) = tmp;
        }
        index_matrices(levels(i)-1) = cfd_idx;
        confd_counts(levels(i)-1) = trans(sum(cfd_idx));
    }

    // move confounding matrices from Rcpp::List into arma::field
    field<mat> cfd_matrices(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        Rcpp::NumericMatrix temp = cfd_factors[i];
        cfd_matrices(i) = mat(temp.begin(), temp.nrow(), temp.ncol(), false);
        row_factor += index_matrices(i) * cfd_matrices(i);
    }

    // find the indices for training and testing elements 
    train_idx = find(train_indicator);
    test_idx = find(train_indicator == 0);

    // check the fitting with initial values
    predict(row_factor, column_factor, cell_factor, gene_factor, predictions);
    residual = data - predictions;
    evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
    loss = compute_loss(cfd_matrices, column_factor, cell_factor, gene_factor, lambda1, lambda2, alpha, sum_residual, 1);

    while(iter <= max_iter) {

        if(iter % 10 == 0){
            cout << "Iteration " << iter << " ---------------------------------" << endl;
        }
        
        // update all confonding matrices
        gram = column_factor * column_factor.t();
        for(i = 0; i < cfd_num; i++){

            residual += index_matrices(i) * cfd_matrices(i) * column_factor;
            if(tuning == 1){
                tune_cfd_row(residual, train_indicator, cfd_matrices(i), column_factor, cfd_indicators.col(i), gram, lambda1, tuning);
            }else if(tuning == 0){
                fit_cfd_row(residual, cfd_matrices(i), column_factor, index_matrices(i), confd_counts(i), gram, lambda1, tuning);
            }else{
                cout << "Parameter tuning should be either 0 or 1!" << endl;
                exit(1);
            }

            // a trick to reduce computational burden
            if(i != cfd_num - 1){
                residual -= index_matrices(i) * cfd_matrices(i) * column_factor;
            }
        }

        // update columm_factor
        row_factor.zeros();
        for(i = 0; i < cfd_num; i ++) {
            row_factor += index_matrices(i) * cfd_matrices(i);
        }

        residual = data - cell_factor * gene_factor;
        optimize_cfd_col(residual, train_indicator, row_factor, column_factor, lambda1, tuning, 30);
        
        // update the base factor and coding for cell decomposition
        residual = data - row_factor * column_factor;
        optimize_base(residual, train_indicator, cell_factor, gene_factor, tuning);
        if(tuning == 1){
            tune_coding(residual, train_indicator, cell_factor, gene_factor, lambda2, alpha, tuning);
        }else{
            fit_coding(residual, cell_factor, gene_factor, lambda2, alpha, tuning);
        }
        
        residual -=  cell_factor * gene_factor;

        // check fitting every 10 steps
        if(iter % 10 == 0){
            pre_loss = loss;
            evaluate(residual, train_idx, test_idx, sum_residual, train_rmse, test_rmse, tuning, iter, 1);
            loss = compute_loss(cfd_matrices, column_factor, cell_factor, gene_factor, lambda1, lambda2, alpha, sum_residual, 1);

            delta_loss = pre_loss - loss;
            cout << "Delta loss for iter " << iter << ":" << delta_loss << endl;

            // if(delta_loss/1000 <= 1e-5){
            //     decay = 1e-5;
            // }else if(delta_loss/1000 <= 1e-4){
            //     decay = 1e-4;
            // }else if(delta_loss/1000 <= 1e-3){
            //     decay = 1e-3;
            // }else if(delta_loss/1000 <= 1e-2){
            //     decay = 1e-2;
            // }else if(delta_loss/1000 <= 1e-1){
            //     decay = 1e-1;
            // }else{
            //     decay = 1.0;
            // }

            if(delta_loss/pre_loss < global_tol){
                break;
            }
        }
        iter++;
    }

    // put the updated confounding matrices into a R list
    List row_matrices;
    for(i = 0; i < cfd_num; i ++) {
        row_matrices["factor" + std::to_string(i)] = cfd_matrices(i);
    }

    return List::create(Named("train_indicator") = train_indicator,
                        Named("row_matrices") = row_matrices,
                        Named("column_factor") = column_factor, 
                        Named("cell_factor") = cell_factor,
                        Named("gene_factor") = gene_factor,
                        Named("test_rmse") = test_rmse, 
                        Named("loss") = loss);
}

// [[Rcpp::export]]
List batch_optimize(const mat& data, const List& cfd_factors, mat column_factor, const umat& cfd_indicators,
                    mat cell_factor, mat gene_factor, const unsigned int num_batch, const unsigned int predefined_batch, const uvec batch_assignment, 
                    const double lambda1 = 0.1, const double lambda2 = 0.01, const double alpha = 1.0,
                    const double global_tol = 1e-10, const double sub_tol = 1e-5, const unsigned int max_iter = 10000){

    cout.precision(12);
    // time_t tic, toc;
    unsigned int i, iter = 0, k, rank = cell_factor.n_cols, cfd_num = cfd_factors.size();
    double loss = std::numeric_limits<double>::max(); // decay = 1.0;
    double pre_loss, delta_loss, train_rmse, sum_residual = 0.0;
    uvec ids, ord;

    uvec batch_ids = unique(batch_assignment);
    vec trace_ZtZ = zeros(num_batch), trace_RtR = zeros(num_batch);

    mat gram, residual, batch_row_factor, batch_cell_factor; // row_factor = zeros(size(cell_factor));
    mat cfd_XtX = zeros(rank, rank), cell_StS = zeros(size(cfd_XtX));
    mat cfd_XtZ = zeros(rank, column_factor.n_cols), cell_StZ = zeros(size(cfd_XtZ));

    cube XtX = zeros(rank, rank, num_batch), StS = zeros(size(XtX));
    cube XtZ = zeros(rank, column_factor.n_cols, num_batch), StZ = zeros(size(XtZ));

    // check whether the number of the confounding matrices is equal to the number of confounding indicators.
    if(cfd_num != cfd_indicators.n_cols){
        cout << "The number of confounding matrices should be the same as the column number of cfd_indicators." << endl;
        exit(1);
    }

    // assign row indices into batches, whose sizes are determined by num_batch
    field<uvec> batches(num_batch);
    if(predefined_batch == 1){
        for(i = 0; i < size(batch_ids); i++){
            batches(i) = find(batch_assignment == batch_ids(i));
        }
    }else{
        batches = generate_batches(data.n_rows, num_batch);
    }

    // place indices of confounders into arma::field for computational consideration
    field<mat> index_matrices(cfd_num);
    field<mat> cfd_cnts(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        uvec levels = unique(cfd_indicators.col(i));
        mat cfd_idx = zeros(cfd_indicators.n_rows, levels.n_elem);

        for(unsigned int k = 0; k < levels.n_elem; k++) {
            vec tmp = cfd_idx.col(k);
            tmp.elem(find(cfd_indicators.col(i) == levels(k))).ones();
            cfd_idx.col(k) = tmp;
        }
        index_matrices(i) = cfd_idx;

        mat tmp = zeros(levels.n_elem, num_batch);
        for(unsigned int k = 0; k < num_batch; k++){
            tmp.col(k) = trans(sum(cfd_idx.rows(batches(k))));
        }
        cfd_cnts(i) = tmp;
    }

    // move confounding matrices from List into field
    field<mat> cfd_matrices(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        Rcpp::NumericMatrix temp = cfd_factors[i];
        cfd_matrices(i) = mat(temp.begin(), temp.nrow(), temp.ncol(), false);
    }

    while(iter <= max_iter) {

        if(iter % 10 == 0){
            cout << "Iteration " << iter << " ---------------------------------" << endl;
        }
        // sum_residual = 0.0;
        ord = randperm(num_batch);
        for(unsigned int b = 0; b < num_batch; b++){
            k = ord(b);
            ids = batches(k);

            batch_row_factor = zeros(ids.n_elem, rank);
            for(i = 1; i < cfd_num; i ++) {
                batch_row_factor += index_matrices(i).rows(ids) * cfd_matrices(i);
            }
            residual = data.rows(ids) - batch_row_factor * column_factor - cell_factor.rows(ids) * gene_factor;

            // tic = time(0);
            // update all confonding matrices
            trace_RtR(k) = 0.0;
            gram = column_factor * column_factor.t();
            for(i = 0; i < cfd_num; i++){
                if(i != 0){
                    residual += index_matrices(i).rows(ids) * cfd_matrices(i) * column_factor;
                }

                fit_cfd_row(residual, cfd_matrices(i), column_factor, index_matrices(i).rows(ids), cfd_cnts(i).col(k), gram, lambda1);

                residual -= index_matrices(i).rows(ids) * cfd_matrices(i) * column_factor;
                trace_RtR(k) += trace(trans(cfd_matrices(i)) * cfd_matrices(i));
            }
            // toc = time(0);
            // cout << "Computing updates for cfd_matrices took "<< difftime(toc, tic) <<" second(s)."<< endl;

            batch_row_factor = zeros(ids.n_elem, rank);
            for(i = 0; i < cfd_num; i ++) {
                batch_row_factor += index_matrices(i).rows(ids) * cfd_matrices(i);
            }
            // row_factor.rows(ids) = batch_row_factor;

            // remove past information to accelerate convergence.
            residual = data.rows(ids) - cell_factor.rows(ids) * gene_factor;
            if(iter > 0){
                cfd_XtX -= XtX.slice(k);
                cfd_XtZ -= XtZ.slice(k);
            }
            XtX.slice(k) = trans(batch_row_factor) * batch_row_factor;
            XtZ.slice(k) = trans(batch_row_factor) * residual;
            cfd_XtX += XtX.slice(k);
            cfd_XtZ += XtZ.slice(k);

            // update columm_factor
            if(iter == 0){
                column_factor = solve(cfd_XtX + lambda1 * (b+1) * eye(rank, rank), cfd_XtZ, solve_opts::likely_sympd);
            }else{
                column_factor = solve(cfd_XtX + lambda1 * num_batch * eye(rank, rank), cfd_XtZ, solve_opts::likely_sympd);
            }

            // update the base factor and coding for cell decomposition
            residual = data.rows(ids) - batch_row_factor * column_factor;
            batch_cell_factor = cell_factor.rows(ids);
            fit_coding(residual, batch_cell_factor, gene_factor, lambda2, alpha, 0);
            cell_factor.rows(ids) = batch_cell_factor;

            // remove past information to improve convergence.
            if(iter > 0){
                cell_StS -= StS.slice(k);
                cell_StZ -= StZ.slice(k);
            }

            StS.slice(k) = trans(batch_cell_factor) * batch_cell_factor;
            StZ.slice(k) = trans(batch_cell_factor) * residual;
            // trace_ZtZ(k) = trace(trans(residual) * residual);
            trace_ZtZ(k) = pow(norm(residual, "F"), 2);
            cell_StS += StS.slice(k);
            cell_StZ += StZ.slice(k);

            if(iter == 0){
                gene_factor = lagrange_dual(cell_StS/(b+1), cell_StZ/(b+1), 1, 100, 1e-4);
            }else{
                gene_factor = lagrange_dual(cell_StS/num_batch, cell_StZ/num_batch, 1, 100, 1e-4);
            }
            // residual -= batch_cell_factor * gene_factor;
            // sum_residual += accu(sum(square(residual), 1));
        }

        // check fitting every 10 steps
        if(iter % 10 == 0){

            // compute the total sum of squared difference
            sum_residual = trace(trans(gene_factor) * cell_StS * gene_factor) - 2 * trace(trans(gene_factor) * cell_StZ) + sum(trace_ZtZ);

            // calculate RMSE on the train set
            train_rmse = std::sqrt(sum_residual/(data.n_rows * data.n_cols));

            pre_loss = loss;

            // add squared residual
            loss = 0.5 * sum_residual/num_batch;

            // add penalty from confounding representation
            loss +=  0.5 * lambda1 * (sum(trace_RtR)/num_batch + trace(trans(column_factor) * column_factor));

            // add penalty on cell representation
            if(alpha != 1){
                loss += 0.5 * lambda2 * (1 - alpha) * pow(norm(cell_factor, "F"), 2)/num_batch;
            }
            loss += lambda2 * alpha * sum(sum(abs(cell_factor), 1))/num_batch;

            delta_loss = std::abs(pre_loss - loss);

            cout << "Train RMSE for iter " << iter << ":" << train_rmse << endl;
            cout << "Loss for iter " << iter << ":" << loss << endl;
            cout << "Delta loss for iter " << iter << ":" << delta_loss << endl;

            // if(delta_loss/1000 <= 1e-5){
            //     decay = 1e-5;
            // }else if(delta_loss/1000 <= 1e-4){
            //     decay = 1e-4;
            // }else if(delta_loss/1000 <= 1e-3){
            //     decay = 1e-3;
            // }else if(delta_loss/1000 <= 1e-2){
            //     decay = 1e-2;
            // }else if(delta_loss/1000 <= 1e-1){
            //     decay = 1e-1;
            // }else{
            //     decay = 1.0;
            // }

            if(delta_loss/pre_loss < global_tol){
                break;
            }
        }
        iter++;
    }

    residual = data - cell_factor * gene_factor;
    for(i = 1; i < cfd_num; i ++) {
        residual -= index_matrices(i) * cfd_matrices(i) * column_factor;
    }

    for(i = 0; i < cfd_num; i++){
        
        if(i != 0){
            residual += index_matrices(i) * cfd_matrices(i) * column_factor;
        }

        vec cfd_cnt = sum(cfd_cnts(i), 1);
        fit_cfd_row(residual, cfd_matrices(i), column_factor, index_matrices(i), cfd_cnt, gram, lambda1);

        // a trick to reduce computational burden
        if(i != cfd_num - 1){
            residual -= index_matrices(i) * cfd_matrices(i) * column_factor;
        }
    }

    // put the updated confounding matrices into a R list
    List row_matrices;
    for(i = 0; i < cfd_num; i ++) {
        row_matrices["factor" + std::to_string(i)] = cfd_matrices(i);
    }

    return List::create(Named("row_matrices") = row_matrices,
                        Named("column_factor") = column_factor,
                        Named("cell_factor") = cell_factor,
                        Named("gene_factor") = gene_factor,
                        Named("loss") = loss);
}

// [[Rcpp::export]]
List sample_optimize(const mat& data, const List& cfd_factors, mat column_factor, const umat& cfd_indicators,
                    mat cell_factor, mat gene_factor, const unsigned int num_batch, const unsigned int predefined_batch, const uvec batch_assignment, 
                    const double lambda1 = 0.1, const double lambda2 = 0.01, const double alpha = 1.0,
                    const double global_tol = 1e-10, const double sub_tol = 1e-5, const unsigned int max_iter = 10000){

    cout.precision(12);
    // time_t tic, toc;
    unsigned int i, iter = 0, k, rank = cell_factor.n_cols, cfd_num = cfd_factors.size();
    double loss = std::numeric_limits<double>::max(); // decay = 1.0;
    double pre_loss, delta_loss, train_rmse, sum_residual = 0.0;
    uvec ids, ord;

    uvec batch_ids = unique(batch_assignment);
    vec trace_ZtZ = zeros(num_batch), trace_RtR = zeros(num_batch);

    mat gram, residual, batch_row_factor, batch_cell_factor; // row_factor = zeros(size(cell_factor));
    mat XtX = zeros(rank, rank), Xty = zeros(rank, column_factor.n_cols), inv_XtX = zeros(rank, rank);
    mat cfd_XtX = zeros(rank, rank), cell_StS = zeros(size(cfd_XtX));
    mat cfd_XtZ = zeros(rank, column_factor.n_cols), cell_StZ = zeros(size(cfd_XtZ));

    cube XtX = zeros(rank, rank, num_batch), StS = zeros(size(XtX));
    cube XtZ = zeros(rank, column_factor.n_cols, num_batch), StZ = zeros(size(XtZ));

    // check whether the number of the confounding matrices is equal to the number of confounding indicators.
    if(cfd_num != cfd_indicators.n_cols){
        cout << "The number of confounding matrices should be the same as the column number of cfd_indicators." << endl;
        exit(1);
    }

    // place indices of confounders into arma::field for computational consideration
    field<mat> index_matrices(cfd_num);
    field<vec> confd_counts(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        uvec levels = unique(cfd_indicators.col(i));
        mat cfd_idx = zeros(cfd_indicators.n_rows, levels.n_elem);

        for(unsigned int k = 0; k < levels.n_elem; k++) {
            vec tmp = cfd_idx.col(k);
            tmp.elem(find(cfd_indicators.col(i) == levels(k))).ones();
            cfd_idx.col(k) = tmp;
        }
        index_matrices(i) = cfd_idx;
        confd_counts(i) = trans(sum(cfd_idx));
    }

    // move confounding matrices from Rcpp::List into arma::field
    field<mat> cfd_matrices(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        Rcpp::NumericMatrix temp = cfd_factors[i];
        cfd_matrices(i) = mat(temp.begin(), temp.nrow(), temp.ncol(), false);
        row_factor += index_matrices(i) * cfd_matrices(i);
    }

    // check the fitting with initial values
    residual = data - row_factor * column_factor;

    // compute loss
    loss = sum_residual/2;
    for(unsigned int i = 0; i < cfd_matrices.n_elem; i++){
        loss += lambda1 * pow(norm(cfd_matrices(i), "F"), 2)/2;
    }
    loss += lambda1 * pow(norm(column_factor, "F"), 2)/2;

    cout << "Begin step 1:" << endl;
    while(iter <= max_iter) {

        if(iter % 10 == 0){
            cout << "Iteration " << iter << " ---------------------------------" << endl;
        }
        
        // update all confonding matrices
        gram = column_factor * column_factor.t();
        for(i = 0; i < cfd_num; i++){

            residual += index_matrices(i) * cfd_matrices(i) * column_factor;
            fit_cfd_row(residual, cfd_matrices(i), column_factor, index_matrices(i), confd_counts(i), gram, lambda1);

            // a trick to reduce computational burden
            if(i != cfd_num - 1){
                residual -= index_matrices(i) * cfd_matrices(i) * column_factor;
            }
        }

        // update columm_factor
        row_factor.zeros();
        for(i = 0; i < cfd_num; i ++) {
            row_factor += index_matrices(i) * cfd_matrices(i);
        }
        
        XtX = row_factor.t() * row_factor;
        Xty = row_factor.t() * residual;
        XtX.diag() += lambda; 
        inv_XtX = inv(XtX);

        #if defined(_OPENMP)
            #pragma omp parallel for num_threads(n_cores) schedule(dynamic, 50)
        #endif
        for(unsigned int i = 0; i < residual.n_cols; i++) {
            c_factor.col(i) = inv_XtX * Xty.col(i);
        }
        residual = data - row_factor * column_factor;

        // check fitting every 10 steps
        if(iter % 10 == 0){

            // compute loss
            pre_loss = loss;
            loss = sum_residual/2;
            for(unsigned int i = 0; i < cfd_matrices.n_elem; i++){
                loss += lambda1 * pow(norm(cfd_matrices(i), "F"), 2)/2;
            }
            loss += lambda1 * pow(norm(column_factor, "F"), 2)/2; 

            delta_loss = pre_loss - loss;
            cout << "Delta loss for iter " << iter << ":" << delta_loss << endl;

            if(delta_loss/pre_loss < global_tol){
                break;
            }
        }
        iter++;
    }

    cout << "Begin step 2:" << endl;

    // assign row indices into batches, whose sizes are determined by num_batch
    field<uvec> batches(num_batch);
    if(predefined_batch == 1){
        for(i = 0; i < size(batch_ids); i++){
            batches(i) = find(batch_assignment == batch_ids(i));
        }
    }else{
        batches = generate_batches(data.n_rows, num_batch);
    }

    // place indices of confounders into arma::field for computational consideration
    field<mat> index_matrices(cfd_num);
    field<mat> cfd_cnts(cfd_num);
    for(i = 0; i < cfd_num; i ++) {
        uvec levels = unique(cfd_indicators.col(i));
        mat cfd_idx = zeros(cfd_indicators.n_rows, levels.n_elem);

        for(unsigned int k = 0; k < levels.n_elem; k++) {
            vec tmp = cfd_idx.col(k);
            tmp.elem(find(cfd_indicators.col(i) == levels(k))).ones();
            cfd_idx.col(k) = tmp;
        }
        index_matrices(i) = cfd_idx;

        mat tmp = zeros(levels.n_elem, num_batch);
        for(unsigned int k = 0; k < num_batch; k++){
            tmp.col(k) = trans(sum(cfd_idx.rows(batches(k))));
        }
        cfd_cnts(i) = tmp;
    }

    while(iter <= max_iter) {

        if(iter % 10 == 0){
            cout << "Iteration " << iter << " ---------------------------------" << endl;
        }
        
        // sum_residual = 0.0;
        ord = randperm(num_batch);
        for(unsigned int b = 0; b < num_batch; b++){
            k = ord(b);
            ids = batches(k);

            batch_cell_factor = cell_factor.rows(ids);
            fit_coding(residual.rows(ids), batch_cell_factor, gene_factor, lambda2, alpha, 0);
            cell_factor.rows(ids) = batch_cell_factor;

            // remove past information to improve convergence.
            if(iter > 0){
                cell_StS -= StS.slice(k);
                cell_StZ -= StZ.slice(k);
            }

            StS.slice(k) = trans(batch_cell_factor) * batch_cell_factor;
            StZ.slice(k) = trans(batch_cell_factor) * residual.rows(ids);
            // trace_ZtZ(k) = trace(trans(residual) * residual);
            // trace_ZtZ(k) = pow(norm(residual, "F"), 2);
            trace_ZtZ(k) = pow(norm(residual.rows(ids), "F"), 2);
            cell_StS += StS.slice(k);
            cell_StZ += StZ.slice(k);

            if(iter == 0){
                gene_factor = lagrange_dual(cell_StS/(b+1), cell_StZ/(b+1), 1, 100, 1e-4);
            }else{
                gene_factor = lagrange_dual(cell_StS/num_batch, cell_StZ/num_batch, 1, 100, 1e-4);
            }
            // residual -= batch_cell_factor * gene_factor;
            // sum_residual += accu(sum(square(residual), 1));
        }

        // check fitting every 10 steps
        if(iter % 10 == 0){

            // compute the total sum of squared difference
            // sum_residual = trace(trans(gene_factor) * cell_StS * gene_factor) - 2 * trace(trans(gene_factor) * cell_StZ) + sum(trace_ZtZ);
            sum_residual = pow(norm(residual - cell_factor * gene_factor, "F"), 2);

            // calculate RMSE on the train set
            train_rmse = std::sqrt(sum_residual/(data.n_rows * data.n_cols));

            pre_loss = loss;

            // add squared residual
            loss = 0.5 * sum_residual/num_batch;

            // add penalty on cell representation
            if(alpha != 1){
                loss += 0.5 * lambda2 * (1 - alpha) * pow(norm(cell_factor, "F"), 2)/num_batch;
            }
            loss += lambda2 * alpha * sum(sum(abs(cell_factor), 1))/num_batch;

            delta_loss = std::abs(pre_loss - loss);

            cout << "Train RMSE for iter " << iter << ":" << train_rmse << endl;
            cout << "Loss for iter " << iter << ":" << loss << endl;
            cout << "Delta loss for iter " << iter << ":" << delta_loss << endl;

            if(delta_loss/pre_loss < global_tol){
                break;
            }
        }
        iter++;
    }

    // put the updated confounding matrices into a R list
    List row_matrices;
    for(i = 0; i < cfd_num; i ++) {
        row_matrices["factor" + std::to_string(i)] = cfd_matrices(i);
    }

    return List::create(Named("row_matrices") = row_matrices,
                        Named("column_factor") = column_factor,
                        Named("cell_factor") = cell_factor,
                        Named("gene_factor") = gene_factor,
                        Named("loss") = loss);
}

