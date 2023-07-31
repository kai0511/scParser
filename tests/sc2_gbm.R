require(SC2)

lambda1 <- 100  # original value 960.5335 
lambda2 <- 0.4 # original value 10, 0.8, 0.5  

# ratio <- 0.1 # ratio for randomly sampling a proportion of observations for model selection
ratio <- 1

setwd("~/data/single_cell_datasets/SCP1985/")
cfd <- read.csv("~/data/single_cell_datasets/SCP1985/gbm_cfd_meta.csv", row.names = 1)
# load("~/data/single_cell_datasets/SCP1985/GBM_2000hvf.RData")
load("~/data/single_cell_datasets/SCP1985/GBM_2000hvf_GSMID.RData")
clusters <- read.csv("SC2_gbm_result.csv", row.names = 1, stringsAsFactors = F, header = T)
clusters <- subset(clusters, louvain %in% c(0, 1, 2, 3, 6))
clusters$louvain[clusters$louvain == 6] <- 4  # rename the 6th cluster

if(ratio <  1){
    data_list <- lapply(data_list, function(d){
        obs_num <- ncol(d)
        idx <- sample(obs_num, floor(obs_num * ratio))
        d[, idx]
    })
}

dataset <- as.matrix(t(do.call(cbind, data_list)))
confounders <- as.matrix(cbind(cfd[rownames(clusters), 1:2], clusters[, 7, drop = F]))
n_cluster <- length(unique(confounders[,3]))

cluster_id <- (confounders[,1] - 1) * n_cluster + 1 + confounders[,3]
confounders <- cbind(confounders, cluster_id)

# sc_factor_num <- 37 
sc_factor_num <- 15 
# sc_factor_num <- seq(5, 30, by = 2)

# object <- SC2(dataset, confounders, split_ratio = 0.1, tuning_iter = 20)
# tune(object, as.integer(sc_factor_num), lambda1 = c(0.1, 1, 5, 10, 20, 50, 100), lambda2 = seq(0.3, 0.7, by = 0.1))
object <- SC2(dataset[rownames(clusters),], confounders[,4, drop = F], split_ratio = 0.1, global_tol = 1e-8, max_iter = 500)
rm(dataset, confounders); gc();
# object <- fit(object, as.integer(sc_factor_num), batch_num = 40, is_batch = T, lambda1 = lambda1, lambda2 = lambda2)
# save(object, file = paste0("SC2_gbm_R", sc_factor_num, "_batch_fitted_object.RData"))

object <- partial_tune(object, cfd_rank = as.integer(sc_factor_num), lambda1 = c(1, 10, 20, 30, 40, 50))
# object <- partial_fit(object, cfd_rank = as.integer(cfd_factor_num), lambda1 = 1)
# save(object, file = paste0("SC2_gbm_R", sc_factor_num, "_fitted_cluster.RData"))
