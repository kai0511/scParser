require(SC2)

lambda1 <- 15.5 # original value 960.5335 
lambda2 <- 0.6 # original value 10, 0.8, 0.5  

# ratio <- 0.05 # ratio for randomly sampling a proportion of observations for model selection
ratio <- 1 

cfd <- read.csv("/home/kai/data/single_cell_datasets/PMC7612735/global_cfd_meta.csv", row.names = 1)
load("~/data/single_cell_datasets/PMC7612735/global_2000hvg.RData")

if(ratio <  1){
    data_list <- lapply(data_list, function(d){
        obs_num <- ncol(d)
        idx <- sample(obs_num, floor(obs_num * ratio))
        d[, idx]
    })
}

dataset <- as.matrix(t(do.call(cbind, data_list)))
confounders <- as.matrix(cfd[rownames(dataset), 4:3])


sc_factor_num <- 23
# sc_factor_num <- seq(15, 40, by = 2)

# object <- SC2(dataset, confounders, split_ratio = 0.1, tuning_iter = 20)
# tune(object, as.integer(sc_factor_num), lambda1 = seq(2, 20, length.out = 5), lambda2 = seq(0.4, 0.8, by = 0.1))
object <- SC2(dataset, confounders, split_ratio = 0.1, tuning_iter = 20, global_tol = 1e-8, max_iter = 500)
rm(dataset, confounders); gc();
object <- fit(object, as.integer(sc_factor_num), batch_num = 60, is_batch = T, lambda1 = lambda1, lambda2 = lambda2)
save(object, file = paste0("SC2_immune_R", sc_factor_num, "_batch_fitted_object.RData"))
