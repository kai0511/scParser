## SR2

SR2 is an ensemble computational tool for interpretive single-cell RNA-seq data analysis. It decomposes variation from multiple biological conditions and cellular variation across bio-samples into shared low-rank latent spaces. Our approach enables various downstream analysis. See the preprint in References section for details of downstream analysis.

## Dependencies
The R package of SR2 has the following dependencies. 
```{r}
Rcpp (>= 1.0.9), RcppArmadillo, dplyr
```

## Installation
SR2 employs **openmp** to support parallel computing. One can install **openmp** before the installation of the package to enjoy the benefit. This package can be installed via the following two ways.

### Install via devtools
```{r}
install.packages("devtools")
install_github("kai0511/SR2")
```
### Install locally
Download the zip file from GitHub, unzip it, and go into the directory. Then, it can be installed by running the following commands
```{Shell}
R CMD build .
R CMD INSTALL SR2_1.0.tar.gz 
```

## Usage

* Data preparation

Here we use the DM dataset (22753*2000) as data matrix as toy example for illustration. The data is place in the data directory of this package.
```{r}
require(SR2)

load("./data/E-HCAD-31_log_transformed_matrix.rds")  # load the data and the exact path for the example data matrix depends.

head(dataset[,1:4])
                        disease_id donnor_id ENSG00000108849 ENSG00000157005
SRR5818088-AAAAAAAAAAAA          1         1        0.000000        0.000000
SRR5818088-AAAAAAAATAGG          1         1        0.000000        2.954825
SRR5818088-AAAAGACGAACG          1         1        7.223458        0.000000
SRR5818088-AAAAGGGCGAAC          1         1        0.000000        2.279514
SRR5818088-AAAAGGGGAACA          1         1        0.000000        0.000000
SRR5818088-AAAAGGGGAACC          1         1        0.000000        0.000000

end_idx <- 2  # The end index for covariate matrix
data[is.na(dataset)] <- 0 # cast NAs to zeros

# In the example data, there are 2 biological variables: diabetes diagnosis (disease_id), donnor id.
confounders <- as.matrix(dataset[ ,1:end_idx])   # matrix for biological variables
```

* Create SR2 object
```{r}
object <- SR2(as.matrix(dataset[,-c(1,2)]), confounders, split_ratio = 0.1, global_tol = 1e-8, sub_tol = 1e-5, tuning_iter = 20)
```
It needs the following arguments:
1. *data*: A log-transformed expression data matrix;
2. *confounder*: A confounder matrix. The elements of the matrix are used as indices to extract corresponding latent representation, so its elements are integer and greater than 0;
3. *split_ratio*: define the proportion of elements in the data matrix used as test set for hyperparameter tuning.  
4. *global_tol*: defines global convergence tolerance for SR2. Note SR2 check convergence every 10 iterations, and the default value for global_tol is 1e-8.
5. *sub_tol*: defines the convergence criteria for elastic net problems. Its impact on global convergence rate is small. By default, it is 1e-5.
6. *tuning_iter*: the number of iterations to run for each try of hyperparameter combinations. In practice, 20 or 30 iterations per try work fine in practice.
7. *max_iter*: the maximum number of iterations. When it reaches, iteration will terminate even if the global convergence criteria do not meet. Its default value is 50000.

* Tune hyperparameters
```{r}
object <- tune(object, latent_dimension = as.integer(seq(10, 30, by = 2)), lambda = seq(1, 20, by = 2), alpha = c(0.2, 0.3, 0.4, 0.5))
```
It needs the following arguments:
1. *object*: An SR2 object created with the above arguments;
2. *latent_dimension*: A integer vector from which the rank of latent dimension is chosen;
3. *lambda*: A numeric vector from which the tuning parameter *lambda* is selected;
4. *alpha*: A numeric vector from which the tuning parameter *alpha* is selected;

* After parameter tuning, the results for tuning will be saved in the current directory. One chose the combination of hyperparameters with the lowest RMSE on test, and fit SR2 with it.
```{r}
# selected hyperparameters for SR2
num_factors <- 23
lambda <- 10
alpha <- 0.4
object <- fit(object, as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = paste0("SR2_ageing_R", num_factors, "_fitted_object.RData"))

> str(object)
List of 7
 $ data           : num [1:377, 1:5000] 0 0 0 0 0.336 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:5000] "X499304660" "X499304661" "X499304664" "X499304666" ...
 $ confounder     : num [1:377, 1:4] 1 1 1 1 1 1 2 2 2 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:4] "pid" "sid" "did" "interaction_indicator"
 $ train_indicator: int [1:377, 1:5000] 1 1 1 1 0 1 1 1 1 1 ...
 $ params         :List of 4
  ..$ global_tol : num 1e-10
  ..$ sub_tol    : num 1e-05
  ..$ tuning_iter: num 30
  ..$ max_iter   : num 50000
 $ cfd_matrices   :List of 4
  ..$ factor0: num [1:2, 1:23] -0.155 0.128 0.144 -0.217 0.862 ...
  ..$ factor1: num [1:8, 1:23] -0.0301 -0.057 -0.2099 0.1279 0.0856 ...
  ..$ factor2: num [1:107, 1:23] -0.377 -0.778 -0.11 -0.552 0.251 ...
  ..$ factor3: num [1:16, 1:23] -0.000318 0.196197 -0.040399 0.031864 0.053713 ...
 $ column_factor  : num [1:23, 1:5000] 0 0.00386 0 -0.00831 0 ...
 - attr(*, "class")= chr "SR2"
```

The fitted object obtained from the above command is an R list object, containing the following elements:
1. log-transformed expression data matrix;
2. confounder matrix
3. train_indicator: an indicator matrix for elements to be concluded as train set.
4. params: parameter setting for SR2
6. cfd_matrices: a list of low-rank representations for biological variables and interaction. One can access the low-rank representation for a specific biological variable with the index of the variable in the confounder matrix. The low-rank representation for the interaction is the last matrix of the cfd_matrices List.
7. column_factor: gene latent representation matrix of K * M, where K is the num_factors and M is the number of genes.

For downstream analysis with results from SR2, please refer to preprint in references. 

## References
Zhao, Kai, et al. "SR2: Interpretable Sparse Matrix Decomposition for Bulk RNA Expression Data Analysis." bioRxiv (2022): 2022-11.

