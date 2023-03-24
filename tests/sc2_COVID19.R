require(SC2)

lambda1 <- 32.86
lambda2 <- 0.4   

# ratio <- 0.05 # ratio for randomly sampling a proportion of observations for model selection
ratio <- 1

load("~/data/single_cell_datasets/single-cell-with-phenotype/E-MTAB-9221-COVID-19/E-MTAB-9221_2000hv.RData")
# load("~/data/single_cell_datasets/single-cell-with-phenotype/E-GEOD-150728-COVID-19/E-GEOD-150728_lognorm_count.RData")

if(ratio < 1){
    donor_id <- unique(dataset[ ,1])
    
    sampled_id <- lapply(donor_id, function(i) {
        row_ids <- which(dataset[ ,1] == i)
        sample(row_ids, floor(length(row_ids) * ratio))
    })
    sampled_id <- do.call(c, sampled_id)
    tuning_data <- dataset[sampled_id, ]
}else{
    tuning_data <- dataset
}

dataset <- as.matrix(tuning_data[,-c(1,2)])
confounders <- cbind(status = tuning_data[,2], donnor_id = tuning_data[,1])

# sc_factor_num <- 12
sc_factor_num <- 18
# sc_factor_num <- seq(10, 30, by = 2)

object <- SC2(dataset, confounders, split_ratio = 0.2, tuning_iter = 20, global_tol = 1e-8)
# tune(object, as.integer(sc_factor_num), lambda1 = seq(20, 50, length.out = 8), lambda2 = c(0.3, 0.4, 0.5, 0.6))
object <- fit(object, as.integer(sc_factor_num), lambda1 = lambda1, lambda2 = lambda2)
save(object, file = paste0("SC2_covid_R", sc_factor_num, "_fitted_object.RData"))
