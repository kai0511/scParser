#' Create an An object with the provided parameters
#'
#' @param object An SC2 oject
#' @param split_ratio A proportion of elements will be drawn from data as testset, and the left over is trainset.
#' @param global_tol The global convergence criteria. One meet, the optimization will be terminated.
#' @param sub_tol The convergence creteria for elastic net problems, which decays with decreasing the loss of objective.
#' @param tuning_iter The number of iterations for each combination of tuning parameters
#' @param max_iter The maxiumme number of iterations. When it is reached, iteraction will terminated even if the global convergence criteria does not meet.
#'
#' @return An SC2 object
#' @export
#' 
SR2 <- function(data, confounders, split_ratio = 0.1, global_tol = 1e-8, sub_tol = 1e-5, tuning_iter = 30, max_iter = 10000){

    # create insider class
    object <- structure(list(), class = "SR2")

    object[['data']] <- data
    object[['confounders']] <- confounders
    object[['split_ratio']] <- split_ratio
    params <- list(global_tol = global_tol, sub_tol = sub_tol, tuning_iter = tuning_iter, max_iter = max_iter)
    object[['params']] <- params
    return(object)
}

#' @export
partial_tune <- function(object, cfd_rank, lambda1 = 0.1){

    # split data into two pieces and generate indicator for easy operation in C++
    train_indicator <- ratio_splitter(object$data, ratio = object[['split_ratio']])

    if(!is.integer(cfd_rank)){
        stop("Error: cfd_rank should be vectors of integer.")
    }

    if(!is.numeric(lambda1)){
        stop("Error: lambda1 should be numeric scalar.")
    }

    param_grid <- cfd_rank
    global_tol <- object[['params']][['global_tol']]
    tuning_iter <- object[['params']][['tuning_iter']]
    rank_tuning <- NULL; reg_tuning <- NULL

    # tune the rank of latent representations
    if(length(param_grid) > 1){
        
        for(i in 1:length(param_grid)){
            cfd_rank = param_grid[i]

            cat('Rank:', param_grid[i], "---------------------------------\n")

            # create a matrix for each covariate and put them into a List
            confounder_list <- lapply(1:ncol(object$confounders), function(i){
                factor_num <- length(unique(object$confounders[,i]))
                matrix(init_parameters(factor_num * cfd_rank), ncol = cfd_rank)
            })
            
            column_factor <- matrix(init_parameters(cfd_rank * ncol(object$data)), nrow = cfd_rank)

            fitted_obj <- partial_optimize(object$data, train_indicator, confounder_list, column_factor, object$confounders, 0.1, 1, global_tol, tuning_iter)

            if(is.null(rank_tuning)){
                rank_tuning <- c(cfd_rank, fitted_obj$test_rmse)
                rank_tuning <- t(as.matrix(rank_tuning))
            }else{
                rank_tuning <- rbind(rank_tuning, c(cfd_rank, fitted_obj$test_rmse))
            }
            write.csv(rank_tuning, file = paste0('SC2_rank_tuning_result.csv'))
        }
    }

    # select the latent rank with the lowest optimal rmse
    if(length(param_grid) > 1){
        idx <- which.min(rank_tuning[, 2])
        cfd_rank <- param_grid[idx, 1]
    }else{
        cfd_rank <- cfd_rank
    }
    
    # tune lambda1 and lambda2 with the selected latent rank
    param_grid <- lambda1
    if(length(param_grid) > 1){

        for(i in seq(length(param_grid))){
            
            cat('parameter grid:', round(param_grid[i], 2), "---------------------------------\n")
            
            lambda1 <- round(param_grid[i], 2)
            
            confounder_list <- lapply(1:ncol(object$confounders), function(i){
                factor_num <- length(unique(object$confounders[,i]))
                matrix(init_parameters(factor_num * cfd_rank), ncol = cfd_rank)
            })
            column_factor <- matrix(init_parameters(cfd_rank* ncol(object$data)), nrow = cfd_rank)

            fitted_obj <- partial_optimize(object$data, train_indicator, confounder_list, column_factor, object$confounders, lambda1, 1, global_tol, tuning_iter)

            if(is.null(reg_tuning)){
                reg_tuning <- c(round(param_grid[i], 2), fitted_obj$test_rmse)
                reg_tuning <- t(as.matrix(reg_tuning))
            }else{
                reg_tuning <- rbind(reg_tuning, c(round(param_grid[i], 2), fitted_obj$test_rmse))
            }
            write.csv(reg_tuning, file = paste0('SC2_R', cfd_rank, '_reg_tuning_result.csv'))
        }
    }
    return(list(rank_tuning = rank_tuning, reg_tuning = reg_tuning))
}

#' @export
partial_fit <- function(object, cfd_rank, lambda1){

    global_tol <- object[['params']][['global_tol']]
    max_iter <- object[['params']][['max_iter']]
    
    confounder_list <- lapply(1:ncol(object$confounder), function(i){
        factor_num <- length(unique(object$confounder[,i]))
        matrix(init_parameters(factor_num * cfd_rank), ncol = cfd_rank)
    })

    column_factor <- matrix(init_parameters(cfd_rank * ncol(object$data)), nrow = cfd_rank)
    train_indicator <- matrix(as.integer(0), nrow = cfd_rank, ncol = cfd_rank);
    
    fitted_obj <- partial_optimize(object$data, train_indicator, confounder_list, column_factor, 
                                   object$confounder, lambda1, 0, global_tol, max_iter)

    object[['cfd_matrices']] <- fitted_obj[['row_matrices']]
    object[['column_factor']] <- fitted_obj[['column_factor']]
    object[['test_rmse']] <- fitted_obj[['test_rmse']]
    return(object)
}

#' Tuning hyperparameters with the provided parameters.
#'
#' @param object An SC2 oject.
#' @param data
#' @param confounder A matrix of dummy variables, which indicates the belonging of each observation for the variables. Here only discreate variables are allowed, such as disease status, gender, tissue, or other biological variables.
#' @param latent_rank  An integer vector, from which the latent rank for cell representations is selected. 
#' @param cfd_rank An integer vector, from which the latent rank for confound representations is selected. 
#' @param lambda1 A numeric vector of L2 penalty for confound representations.
#' @param lambda2 A numeric vector of L2 penalty for cell representations. Here lambda2 and alpha form parameters for elastic net penalty.
#' @param alpha A double for L1 penalty for cell representations. By default, alpha = 1. In practice, alpha is opt to 1 when tuning with grid search in single-cell analysis. Thus, alpha is a number of double ranging from 0 to 1 and will not be tuned.
#' @return A list of matrices of tuning results.
#' @export
#' 
tune <- function(object, latent_rank, cfd_rank = NULL, lambda1 = 0.1, lambda2 = 0.1, alpha = 1){

    # split data into two pieces and generate indicator for easy operation in C++
    train_indicator <- ratio_splitter(object$data, ratio = object[['split_ratio']])

    if((!is.integer(cfd_rank) & !is.null(cfd_rank)) | !is.integer(latent_rank)){
        stop("Error: cfd_rank and latent_rank should be vectors of integer. cfd_rank can be NULL when unspecified.")
    }

    if(!is.numeric(lambda1) | !is.numeric(lambda2) | !is.numeric(alpha)){
        stop("Error: lambda1, lambda2, and alpha should be (vectors of) numeric.")
    }

    if(is.null(cfd_rank)){
        param_grid <- data.frame(cfd_rank = latent_rank, latent_rank = latent_rank)
    }else{
        param_grid <- expand.grid(cfd_rank = cfd_rank, latent_rank = latent_rank)
    }

    global_tol <- object[['params']][['global_tol']]
    sub_tol <- object[['params']][['sub_tol']]
    tuning_iter <- object[['params']][['tuning_iter']]
    rank_tuning <- NULL; reg_tuning <- NULL

    # tune the rank of latent representations
    if(nrow(param_grid) > 1){
        
        for(i in 1:nrow(param_grid)){
            cfd_rank = param_grid[i, 1]
            latent_rank = param_grid[i, 2]

            cat('Rank:', paste(param_grid[i,], collapse = ","), "---------------------------------\n")

            # create a matrix for each covariate and put them into a List
            confounder_list <- lapply(1:ncol(object$confounders), function(i){
                factor_num <- length(unique(object$confounders[,i]))
                matrix(init_parameters(factor_num * cfd_rank), ncol = cfd_rank)
            })
            
            column_factor <- matrix(init_parameters(cfd_rank * ncol(object$data)), nrow = cfd_rank)
            cell_factor <- matrix(init_parameters(latent_rank * nrow(object$data)), ncol = latent_rank)
            gene_factor <- matrix(init_parameters(latent_rank * ncol(object$data)), nrow = latent_rank)
            
            # rescale each row of gene_factor to have unit l2 norm
            row_std <- apply(gene_factor, 1, function(x) sqrt(sum(x^2)))
            gene_factor <- gene_factor/row_std

            fitted_obj <- optimize(object$data, train_indicator, confounder_list, column_factor, object$confounders, cell_factor, gene_factor, 0.1, 0.01, alpha, 1, global_tol, sub_tol, tuning_iter)

            if(is.null(rank_tuning)){
                rank_tuning <- c(cfd_rank, latent_rank, fitted_obj$test_rmse)
                rank_tuning <- t(as.matrix(rank_tuning))
            }else{
                rank_tuning <- rbind(rank_tuning, c(cfd_rank, latent_rank, fitted_obj$test_rmse))
            }
            write.csv(rank_tuning, file = paste0('SC2_rank_tuning_result.csv'))
        }
    }

    # select the latent rank with the lowest optimal rmse
    if(nrow(param_grid) > 1){
        idx <- which.min(rank_tuning[ ,3])
        cfd_rank <- param_grid[idx, 1]
        latent_rank <- param_grid[idx, 2]
    }else{
        cfd_rank <- param_grid[1,1]
        latent_rank <- param_grid[1,2]
    }
    
    # tune lambda1 and lambda2 with the selected latent rank
    param_grid <- expand.grid(lambda1 = lambda1, lambda2 = lambda2)
    if(nrow(param_grid) > 1){

        for(i in seq(nrow(param_grid))){
            
            cat('parameter grid:', paste(round(param_grid[i,], 2), collapse = ','), "---------------------------------\n")
            
            lambda1 <- round(param_grid[i, 1], 2)
            lambda2 <- round(param_grid[i, 2], 2)
            
            confounder_list <- lapply(1:ncol(object$confounders), function(i){
                factor_num <- length(unique(object$confounders[,i]))
                matrix(init_parameters(factor_num * cfd_rank), ncol = cfd_rank)
            })

            column_factor <- matrix(init_parameters(cfd_rank* ncol(object$data)), nrow = cfd_rank)
            cell_factor <- matrix(init_parameters(latent_rank * nrow(object$data)), ncol = latent_rank)
            gene_factor <- matrix(init_parameters(latent_rank * ncol(object$data)), nrow = latent_rank)

            fitted_obj <- optimize(object$data, train_indicator, confounder_list, column_factor, object$confounders,
                                   cell_factor, gene_factor, lambda1, lambda2, alpha, 1, global_tol, sub_tol, tuning_iter)

            if(is.null(reg_tuning)){
                reg_tuning <- c(round(param_grid[i, ], 2), fitted_obj$test_rmse)
                reg_tuning <- t(as.matrix(reg_tuning))
            }else{
                reg_tuning <- rbind(reg_tuning, c(round(param_grid[i,], 2), fitted_obj$test_rmse))
            }
            write.csv(reg_tuning, file = paste0('SC2_R', latent_rank, '_reg_tuning_result.csv'))
        }
    }
    return(list(rank_tuning = rank_tuning, reg_tuning = reg_tuning))
}

#' @export
fit <- function(object, latent_rank, cfd_rank = NULL, batch_num = NULL, is_batch = FALSE, batch_assignment = NULL, lambda1 = 1, lambda2 = 0.5, alpha = 1){
   
    if(is.null(cfd_rank)){
        cfd_rank <- latent_rank
    }

    if(is_batch & is.null(batch_num)){
        stop("If is_batch is TURE then batch_num should be provided.")
    } 

    global_tol <- object[['params']][['global_tol']]
    sub_tol <- object[['params']][['sub_tol']]
    max_iter <- object[['params']][['max_iter']]
    
    confounder_list <- lapply(1:ncol(object$confounder), function(i){
        factor_num <- length(unique(object$confounder[,i]))
        matrix(init_parameters(factor_num * cfd_rank), ncol = cfd_rank)
    })

    column_factor <- matrix(init_parameters(cfd_rank * ncol(object$data)), nrow = cfd_rank)
    cell_factor <- matrix(init_parameters(latent_rank * nrow(object$data)), ncol = latent_rank)
    gene_factor <- matrix(init_parameters(latent_rank * ncol(object$data)), nrow = latent_rank)

    if(!is_batch){
        train_indicator <- matrix(as.integer(0), nrow = latent_rank, ncol = latent_rank);
        fitted_obj <- optimize(object$data, train_indicator, confounder_list, column_factor, object$confounder,
                               cell_factor, gene_factor, lambda1, lambda2, alpha, 0, global_tol, sub_tol, max_iter)
    }else{

        if(is.null(batch_assignment)){
            fitted_obj <- batch_optimize(object$data, confounder_list, column_factor, object$confounder,
                                     cell_factor, gene_factor, batch_num, 0, rep(0, nrow(object$data)), lambda1, lambda2, alpha, global_tol, sub_tol, max_iter)
        } else {
            fitted_obj <- batch_optimize(object$data, confounder_list, column_factor, object$confounder,
                                     cell_factor, gene_factor, batch_num, 1, batch_assignment, lambda1, lambda2, alpha, global_tol, sub_tol, max_iter)
        }
    }

    object[['cfd_matrices']] <- fitted_obj[['row_matrices']]
    object[['column_factor']] <- fitted_obj[['column_factor']]
    object[['cell_factor']] <- fitted_obj[['cell_factor']]
    object[['gene_factor']] <- fitted_obj[['gene_factor']]
    object[['test_rmse']] <- fitted_obj[['test_rmse']]
    return(object)
}
