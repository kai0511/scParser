require(SC2)

lambda1 <- 10 
lambda2 <- 0.5   

# ratio <- 0.05 # ratio for randomly sampling a proportion of observations for model selection
ratio <- 1

load("~/data/single_cell_datasets/single-cell-with-phenotype/E-CURD-114_smoking/E-CURD-114_2000hv.RData")

if(ratio < 1){
    donor_id <- unique(data[ ,1])
    
    sampled_id <- lapply(donor_id, function(i) {
        row_ids <- which(dataset[ ,1] == i)
        sample(row_ids, floor(length(row_ids) * ratio))
    })
    sampled_id <- do.call(c, sampled_id)
    tuning_data <- data[sampled_id, ]
}else{
    tuning_data <- data
}

dataset <- as.matrix(tuning_data[,-c(1,2)])
confounders <- cbind(status = tuning_data[,2], donnor_id = tuning_data[,1])

# sc_factor_num <- 12
sc_factor_num <- 22
# sc_factor_num <- seq(15, 25, by = 1)

object <- SC2(dataset, confounders, split_ratio = 0.1, tuning_iter = 20, global_tol = 1e-8)
# tune(object, as.integer(sc_factor_num), lambda1 = c(10, 20, 30, 40), lambda2 = c(0.4, 0.5, 0.6, 0.7))
object <- fit(object, as.integer(sc_factor_num), lambda1 = lambda1, lambda2 = lambda2)
save(object, file = paste0("SC2_smoking_R", sc_factor_num, "_fitted_object.RData"))
