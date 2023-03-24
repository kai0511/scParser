require(SC2)

lambda1 <- 40   # original value 960.5335 
lambda2 <- 0.7  # original value 10, 0.8, 0.5  


# load("~/data/single_cell_datasets/single-cell-with-phenotype/E-HCAD-31-fang2019single-DM/E-HCAD-31_log_transformed_2000features.rds")
# setwd('~/data/single_cell/E-HCAD-31-DM/')
setwd("~/data/single_cell_datasets/single-cell-with-phenotype/E-HCAD-31-fang2019single-DM")
load("E-HCAD-31_log_transformed_matrix.rds")
# load("E-HCAD-31_log_transformed_2000features.rds")
clusters <- read.csv("louvain_cluster_result.csv", header = T, row.names = 1)

d <- cbind(dataset[, 1, drop = F], clusters[rownames(dataset), , drop = F])
n_cluster <- length(unique(d[,2]))

cluster_id <- sapply(1:nrow(d), function(i) if(d[i,1] == 1) {d[i,2]+1} else{(d[i,1] - 1) * n_cluster + 1 + d[i,2]})
confounders <- as.matrix(cluster_id, ncol = 1)

# confounders <- cbind(disease_id = data[,1], donnor_id = data[,2])

cfd_factor_num <- 11 # number of latent factor for DM data
sc_factor_num <- 11 

object <- SC2(as.matrix(dataset[,-c(1,2)]), confounders, split_ratio = 0.1, tuning_iter = 20)
# tune(object, as.integer(cfd_factor_num), as.integer(sc_factor_num), seq(0.01, 1, length.out = 5), seq(0.1, 0.9, by = 0.1))
# tune(object, latent_rank = as.integer(sc_factor_num), lambda1 = c(1, 10, 20, 30, 40, 50), lambda2 = c(0.5, 0.6, 0.7, 0.8, 0.9))
# object <- fit(object, latent_rank = as.integer(sc_factor_num), lambda1 = lambda1, lambda2 = lambda2)
# save(object, file = paste0("SC2_dm_R", sc_factor_num, "_fitted_integrated.RData"))

object <- partial_tune(object, cfd_rank = as.integer(cfd_factor_num), lambda1 = c(1, 10, 20, 30, 40, 50))
# object <- partial_fit(object, cfd_rank = as.integer(cfd_factor_num), lambda1 = 1)
# save(object, file = paste0("SC2_dm_R", sc_factor_num, "_fitted_cluster.RData"))
