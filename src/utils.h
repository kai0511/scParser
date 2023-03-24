#ifndef __UTILS__
#define __UTILS__

#include "../inst/include/SR2_types.h"

field<uvec> generate_batches(const unsigned int& sample_size, const unsigned int& num_batch);

double objective(const mat& X, const vec& y, const vec& beta,
                 const double& lambda, const double& alpha);

double compute_loss(const vec& residual, const vec& beta, const double& lambda, const double& alpha);

void predict(const mat& row_factor, const mat& column_factor, const mat& cell_factor, const mat& gene_factor, mat& predictions);

void evaluate(mat& residual, const uvec& train_idx, const uvec& test_idx, double& sum_residual, 
              double& train_rmse, double& test_rmse, const int& tuning, 
              const int& iter /* = 0*/, const int& verbose /*= 1*/);

double compute_loss(const field<mat>& cfd_factor, const mat& column_factor, const mat& cell_factor, const mat& gene_factor, 
                    const double& lambda1, const double& lambda2, const double& alpha, double& sum_residual, const int& verbose /* = 1 */);

#endif
