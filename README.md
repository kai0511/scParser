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

### Data preparation

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

### Create SR2 object
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
7. *max_iter*: the maximum number of iterations. When it reaches, iteration will terminate even if the global convergence criteria do not meet. Its default value is 10000.

### Tune hyperparameters
```{r}
object <- tune(object, latent_rank = as.integer(seq(10, 30, by = 2)), lambda1 = c(0.1, 1, 10, 20, 50, 100), lambda2 = c(0.01, 0.1, 0.2, 0.4, 0.9))
```
It has the following arguments:
1. *object*: An SR2 object created with the above arguments;
2. *latent_rank*: An integer vector from which the rank of latent dimension is chosen;
3. *cfd_rank* (optional): An integer vector from which the rank of latent dimension for modeling variation from biological conditions is chosen. If not applied, we assume that $K_1 = K_2$. See reference for details.

4. *lambda1*: A numeric vector from which the tuning parameter *lambda1* is selected. $\lambda_1$ controls penalty for latent representations for biological conditions.
5. *lambda2*: A numeric vector from which the tuning parameter *lambda2* is selected. $\lambda_2$ controls penalty for latent representations for cellular representations.
6. *alpha* (optional): A numeric vector from which the tuning parameter *alpha* is selected. By default, $\alpha$ is 1.

* After parameter tuning, the results for tuning will be saved in the current directory. One chose the combination of hyperparameters with the lowest RMSE on test, and fit SR2 with it.

### Model fitting

```{r}
# selected hyperparameters for SR2
num_factors <- 11
lambda1 <- 40   
lambda2 <- 0.7

object <- fit(object, as.integer(num_factors), lambda = lambda, alpha = alpha)
save(object, file = paste0("SR2_DM_R", num_factors, "_fitted_object.RData"))

> str(object)
List of 9
 $ data           : num [1:22753, 1:2000] 0 0 7.22 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:22753] "SRR5818088-AAAAAAAAAAAA" "SRR5818088-AAAAAAAATAGG" "SRR5818088-AAAAGACGAACG" "SRR5818088-AAAAGGGCGAAC" ...
  .. ..$ : chr [1:2000] "ENSG00000108849" "ENSG00000157005" "ENSG00000118785" "ENSG00000164692" ...
 $ train_indicator: int [1:22753, 1:2000] 1 1 0 1 0 1 1 1 1 1 ...
 $ confounder     : num [1:22753, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:22753] "SRR5818088-AAAAAAAAAAAA" "SRR5818088-AAAAAAAATAGG" "SRR5818088-AAAAGACGAACG" "SRR5818088-AAAAGGGCGAAC" ...
  .. ..$ : chr [1:2] "disease_id" "donnor_id"
 $ params         :List of 4
  ..$ global_tol : num 1e-08
  ..$ sub_tol    : num 1e-05
  ..$ tuning_iter: num 30
  ..$ max_iter   : num 10000
 $ cfd_matrices   :List of 2
  ..$ factor0: num [1:2, 1:11] -1.3428 -0.9747 -0.736 -0.7037 0.0449 ...
  ..$ factor1: num [1:9, 1:11] 0.282 0.282 -0.188 -0.11 0.123 ...
 $ column_factor  : num [1:11, 1:2000] -0.1435 -0.15236 -0.13395 0.00664 0.12703 ...
 $ cell_factor    : num [1:22753, 1:11] 0 -4.897 -0.172 0 0.597 ...
 $ gene_factor    : num [1:11, 1:2000] 0.00835 -0.0087 0.01543 0.01838 0.0063 ...
 - attr(*, "class")= chr "SR2"
```

The fitted object obtained from the above command is an R list object, containing the following elements:
1. log-transformed expression data matrix;
2. train_indicator: an indicator matrix for elements to be concluded as train set.
3. confounder matrix
4. params: parameter setting for SR2
6. cfd_matrices: a list of low-rank representations for biological variables. One can access the low-rank representation for a specific biological variable with the index of the variable in the confounder matrix.
7. column_factor: gene latent representation matrix of K * M, where K is the num_factors and M is the number of genes.
8. cell_factor: sparse representation for cells 
9. gene_factor: gene latent representation matrix of K * M, where K is the num_factors and M is the number of genes.

### Modeling the effect of biological conditions on gene expression for cell populations

```{r}
object <- SR2(as.matrix(dataset[,-c(1,2)]), confounders, split_ratio = 0.1, tuning_iter = 20)
```
1. *confounders*: A confounder matrix. The elements of the matrix are used as indices to extract corresponding latent representation, so its elements are integer and greater than 0. To model the interaction between cell populations and biological conditions (e.g., disease status), the confounder matrix should contain a column with each element representing the combination of the cell population and biological condition the corresponding cell belongs. For example, if there are 5 cell populations and 2 levels for a biological condition, then each element of the column should take a value of integer from 1 to 10 to represent the combination it belongs to.
2. Other parameters available is the same as we introduce previously, and the definition of the parameters is also the same.


To tune the parameters for the new model, the following can be employed.
```{r}
object <- partial_tune(object, cfd_rank = as.integer(cfd_factor_num), lambda1 = c(0.1, 1, 10, 20, 30, 40, 50))
```
1. *cfd_rank*: An integer vector from which the rank of latent dimension for modeling variation from biological conditions is chosen.
2. *lambda*: A numeric vector from which the tuning parameter *lambda* is selected. $\lambda$ controls penalty for latent representations for biological variables.

Then, we can fit the new model with the parameters that performs the best in the previous step.
```{r}
object <- partial_fit(object, cfd_rank = as.integer(11), lambda1 = 1)
```

The fitted object obtained from the above command is an R list object as follows:
```{r}
> str(object)
List of 6
 $ data         : num [1:22753, 1:2000] 0 0 7.22 0 0 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:22753] "SRR5818088-AAAAAAAAAAAA" "SRR5818088-AAAAAAAATAGG" "SRR5818088-AAAAGACGAACG" "SRR5818088-AAAAGGGCGAAC" ...
  .. ..$ : chr [1:2000] "ENSG00000108849" "ENSG00000157005" "ENSG00000118785" "ENSG00000164692" ...
  $ confounder     : num [1:22753, 1:2] 3 2 6 5 5 5 5 5 2 4 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:22753] "SRR5818088-AAAAAAAAAAAA" "SRR5818088-AAAAAAAATAGG" "SRR5818088-AAAAGACGAACG" "SRR5818088-AAAAGGGCGAAC" ...
  .. ..$ : chr [1:2] "interaction_id" "donnor_id"
 $ split_ratio  : num 0.1
 $ params       :List of 4
  ..$ global_tol : num 1e-08
  ..$ sub_tol    : num 1e-05
  ..$ tuning_iter: num 20
  ..$ max_iter   : num 50000
 $ cfd_matrices :List of 1
  ..$ factor0: num [1:12, 1:11] -0.537 -0.985 -0.39 0.74 -1.168 ...
  ..$ factor1: num [1:9, 1:11] 0.272 0.241 -0.158 -0.08 0.103 ...
 $ column_factor: num [1:11, 1:2000] 1.527 -0.611 -0.673 1.339 0.426 ...
 - attr(*, "class")= chr "SR2"
```
1. *cfd_matrices*: the first element of the *cfd_matrices* is the latent representations for different combinations of the cell population and biological conditions. 
2. The meaning of other elements is the same as we introduced previously.

For Details of downstream analysis with results from SR2, please refer to preprint in references. 

## References
Zhao, Kai, et al. "SR2: Sparse Representation Learning for Scalable Single-cell RNA Sequencing Data Analysis." bioRxiv (2023): 2023-08.

