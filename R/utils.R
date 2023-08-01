`%+%` <- function(m, v){
    # add vector v to matrix m by column
    if(class(v) == 'matrix'){
        v <- as.vector(v)
    }  
    sweep(m, 2, v, '+')
}

`%-%` <- function(m, v){
    # minus vector v from matrix m by column
    if(class(v) == 'matrix'){
        v <- as.vector(v)
    } 
    sweep(m, 2, v, '-')
}

calculate_idx <- function(idx, num_row){
    # calculate the position in trainset matrix
    col_idx <- floor(idx / num_row)  # column idx
                
    if(idx == col_idx * num_row){
        row_idx <- num_row
    } else {
        row_idx <- idx - col_idx * num_row
        col_idx <- col_idx + 1
    }
    return(c(row_idx, col_idx))
}

init_parameters <- function(size, init_mean = 0.0, init_std = 0.01) {
    return(rnorm(size, mean = init_mean, sd = init_std))
}

split_str <- function(s){
    l <- unlist(strsplit(s, split=c('_')))
    idx <- which(l == 'v7')
    len <- length(l)
    
    disease <- l[1]
    tissue <- paste(l[(idx + 1): len], collapse = '_')
    c(disease, tissue)
}

obtain_indication_matrix <- function(trainset, only_positive = F){
    indication_matrix <- matrix(0, nrow = nrow(trainset), ncol = ncol(trainset))
    if(only_positive){
        indication_matrix[!is.na(trainset)] <- 1
    } else {
        indication_matrix[!is.na(trainset)] <- 1
        indication_matrix[trainset < 0] <- -1
    }
    return(indication_matrix)
}

#' matrix splitting via element-wise sampling without replacement
#'
#' @param data data matrix for splitting
#' @param ratio a proportion of elements considered as testset
#' @param rm.na.col whether remove columns with all zeros
#'
#' @return a list
#' @export
#'
#' @examples results <- ratio_splitter(data)

ratio_splitter <- function(data, ratio = 0.2, rm.na.col = T){
    # default ratio for test is 0.1, that is, 10% obs. will be randomly assigned to testset

    train_indicator <- matrix(as.integer(1), nrow = nrow(data), ncol = ncol(data))
    test_idx <- sample(seq(length(data)), floor(length(data) * ratio), replace = F)
    train_indicator[test_idx] <- as.integer(0)
    return(train_indicator)
}

#' @export 
quantile_norm <- function(cell_factor, clusters, min_cells = 20, quantiles = 50){

    num_clusters <- unique(clusters)
    # use the sample with the largest number of observations as reference dataset
    ref_id <- which.max(unlist(lapply(cell_factor, nrow)))
    
    for (k in seq(length(cell_factor))) {

        for (j in seq(num_clusters)) {

            cells2 <- which(clusters[[k]] == j)
            cells1 <- which(clusters[[ref_id]] == j)
            
            for (i in seq(ncol(cell_factor[[ref_id]]))) {
                num_cells2 <- length(cells2)
                num_cells1 <- length(cells1)
                
                if (num_cells1 < min_cells | num_cells2 < min_cells) {next}
                
                if (num_cells2 == 1) {
                    cell_factor[[k]][cells2, i] <- mean(cell_factor[[ref_id]][cells1, i])
                    next
                }

                q2 <- quantile(cell_factor[[k]][cells2, i], seq(0, 1, by = 1 / quantiles))
                q1 <- quantile(cell_factor[[ref_id]][cells1, i], seq(0, 1, by = 1 / quantiles))
                
                if (sum(q1) == 0 | sum(q2) == 0 | length(unique(q1)) < 2 | length(unique(q2)) < 2) {
                    new_vals <- rep(0, num_cells2)
                } else {
                    warp_func <- stats::approxfun(q2, q1, rule = 2)
                    new_vals <- warp_func(cell_factor[[k]][cells2, i])
                }
                cell_factor[[k]][cells2, i] <- new_vals
            }
        }
    }
    return(cell_factor)
}

#' @export
sc_clustering <- function(adata, cluster_num, method="louvain", min_res=0, max_res=3, max_iter=100){
    
    iter <- 0
    while(iter <= max_iter){

        res <- min_res + (max_res - min_res)/2

        cat("Iter", iter, ":\n")
        
        if(method == "louvain"){
            sc$tl$louvain(adata, resolution = res)
            class_num <- nrow(unique(adata$obs['louvain']))
        }else if(method == "leiden"){
            sc$tl$leiden(adata, resolution = res)
            class_num <- nrow(unique(adata$obs['leiden']))
        }else{
            stop("ERROR: the method parameter only accepts a value from 'louvain' or 'leiden'.")
        }
        cat(class_num, "found at resolution", res, "\n")

        if(class_num > cluster_num){
            max_res <- res
        }else if(class_num < cluster_num){
            min_res <- res
        }else{
            break
        }
        iter <- iter + 1        
    }
    
    cat('INFO:', class_num, 'found at resolution', res, '\n')
    return(list(res, adata))
}

#' @export
mnn_clustering <- function(cell_factor, sample_info, neig_num = 20){
    
    sample_num <- length(sample_info)
    sample_factors <- lapply(sample_info, function(x) cell_factor[x,])
    sample_sizes <- c(0, sapply(sample_factors, nrow))

    edge_list <- lapply(1:sample_num, function(s1){

        edges <- lapply(1:sample_num, function(s2){

            k <- min(neig_num, nrow(sample_factors[[s2]]))
            nns <- RANN::nn2(data = sample_factors[[s2]], 
                             query = sample_factors[[s1]], k, 
                             treetype = "kd", searchtype = "standard")

            start <- rep((sum(sample_sizes[1:s1])+1):sum(sample_sizes[1:(s1+1)]), each = k)
            end <- as.numeric(t(nns$nn.idx)) + sum(sample_sizes[1:s2])
            weights <- as.numeric(t(nns$nn.dists))

            data.frame(from = start, to = end, weight = weights)

        })
        do.call(rbind, edges)
    })

    edges = Reduce(rbind, edge_list)
    ig = igraph::graph_from_data_frame(edges, directed = TRUE)
    ig = igraph::as.undirected(ig, mode = "mutual", edge.attr.comb="first")
    ig = igraph::simplify(ig)
    E(ig)$weight = max(E(ig)$weight) - E(ig)$weight

    cl_louvain = as.factor(igraph::cluster_louvain(ig)$membership)
    clusters = split(cl_louvain, rep(1:length(cell_factor_list), time = nc[-1]))
}


