# Utility functions for scINSIGHT methods. Some published, some not.

#' Normalize \eqn{W_2} and cluster cells.
#' @description
#' Quantile normalization and clustering for the list of \eqn{W_2} (expression matrices of common gene modules). Use weights in Louvain.
#'
#' @param W_2 List of \eqn{W_2}.
#' @param Knn The maximum number of nearest neighbors to search (default 20).
#'
#' @return List of normalized \eqn{W_2}.
#'
#' @importFrom RANN nn2
#' @import igraph
#' @importFrom stats approxfun
#'
#' @examples
#' \dontrun{
#' # Requires scINSIGHT object
#' # Get factorization using K from 5 to 15 and LDA from 0.001 to 10, repeat 5 times
#' # (default setting, can be adjusted for ideal results)
#' scINSIGHTex <- run_scINSIGHT(scINSIGHTex, K = seq(5,15,2),
#'                              LDA = c(0.001, 0.01, 0.1, 1, 10), B = 5)
#' }

norm_clust_strict_weighted = function(W2, Knn = 20, ...){
  ### cosine normalization
  W2norm = lapply(W2, function(x){
    l2 = sqrt(rowSums(x^2))
    l2[l2 < 1e-10] = 1e-10
    x = sweep(x, 1, 1/l2, FUN = "*")
    return(x)
  })
  ### cell number
  nc = sapply(W2, nrow)
  nc = c(0,nc)
  ### find knns
  L = length(W2norm)
  find_knns = lapply(1:L, function(l1){
    tp = lapply(1:L, function(l2){
      nns = RANN::nn2(data = W2norm[[l2]], query = W2norm[[l1]],
                      min(Knn, nrow(W2norm[[l2]])),
                      treetype = "kd", searchtype = "standard")
      from = rep((sum(nc[1:l1])+1):sum(nc[1:(l1+1)]), each = min(Knn, nrow(W2norm[[l2]])))
      to = as.numeric(t(nns$nn.idx)) + sum(nc[1:l2])
      weight = as.numeric(t(nns$nn.dists))
      da = data.frame(from = from, to = to, weight = weight)
      return(da)
    })
    return(Reduce(rbind, tp))
  })
  adjmat = Reduce(rbind, find_knns)
  ig = graph_from_data_frame(adjmat, directed = TRUE)
  # ig = subgraph.edges(ig, eids = E(ig)[which_mutual(ig)], delete.vertices = F)
  ig = as.undirected(ig, mode = "mutual", edge.attr.comb="first")
  ig = simplify(ig)
  E(ig)$weight = max(E(ig)$weight) - E(ig)$weight

  cl_louvain = cluster_louvain(ig)$membership
  clusters = split(cl_louvain, rep(1:L, time = nc[-1]))
  Klv = length(unique(cl_louvain))
  ### normalization
  quantiles = 10
  for(j in 1:Klv){
    ns = sapply(clusters, function(cl) sum(cl==j))
    if(sum(ns > 5) < 2) {next}
    refid = min(which(ns > 5))

    cells1 = which(clusters[[refid]] == j)
    num_cells1 = length(cells1)
    #if (num_cells1 < 2) {next}
    for(l in 1:L){
      if(l == refid) next
      cells2 = which(clusters[[l]] == j)
      num_cells2 = length(cells2)
      if (num_cells2 < 2) {next}
      for(i in 1:ncol(W2[[refid]])){
        q2 = quantile(W2[[l]][cells2, i], seq(0, 1, by = 1 / quantiles))
        q1 = quantile(W2[[refid]][cells1, i], seq(0, 1, by = 1 / quantiles))
        if (sum(q1) == 0 | sum(q2) == 0 | length(unique(q1)) <
            2 | length(unique(q2)) < 2) {
          new_vals = rep(0, num_cells2)
        }
        else {
          warp_func = approxfun(q2, q1, rule = 2)
          new_vals = warp_func(W2[[l]][cells2, i])
        }
        W2[[l]][cells2, i] = new_vals
      }
    }

  }
  return(list(W2=W2, clusters = clusters))
}