require(RANN)


louvainCluster <- function(object, resolution = 1.0, k = 20, prune = 1 / 15, eps = 0.1, nRandomStarts = 10,
                           nIterations = 100, random.seed = 1, verbose = TRUE, dims.use = NULL) {
  output_path <- paste0('edge_', sub('\\s', '_', Sys.time()), '.txt')
  output_path = sub(":","_",output_path)
  output_path = sub(":","_",output_path)
  
  if (is.null(dims.use)) {
    use_these_factors <- 1:ncol(object@H[[1]])
  } else {
    use_these_factors <- dims.use
  }
  
  if (dim(object@H.norm)[1] == 0){
    if (verbose) {
      message("Louvain Clustering on unnormalized cell factor loadings.")
    }
    knn <- RANN::nn2(Reduce(rbind, object@H)[,use_these_factors], k = k, eps = eps)
  } else {
    if (verbose) {
      message("Louvain Clustering on quantile normalized cell factor loadings.")
    }
    knn <- RANN::nn2(object@H.norm[,use_these_factors], k = k, eps = eps)
  }
  snn <- ComputeSNN(knn$nn.idx, prune = prune)
  WriteEdgeFile(snn, output_path, display_progress = FALSE)
  clusts <- RunModularityClusteringCpp(snn,
                                       modularityFunction = 1, resolution = resolution, nRandomStarts = nRandomStarts,
                                       nIterations = nIterations, algorithm = 1, randomSeed = random.seed, printOutput = FALSE,
                                       edgefilename = output_path
  )
  names(clusts) = rownames(object@cell.data)
  rownames(snn) = rownames(object@cell.data)
  colnames(snn) = rownames(object@cell.data)
  clusts <- GroupSingletons(ids = clusts, SNN = snn, verbose = FALSE)
  object@clusters = as.factor(clusts)
  unlink(output_path)
  return(object)
}