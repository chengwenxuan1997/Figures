# # EvaluationMetrics
# 
# if (!grepl("/home/anaconda3/envs/r-4.1.1/bin", Sys.getenv("PATH"))){
#   Sys.setenv("PATH" = paste0(Sys.getenv("PATH"), ":/home/anaconda3/envs/r-4.1.1/bin"))
# }

Rcpp::sourceCpp(file.path(code.path, "EvaluationMetrics.cpp"))
Rcpp::sourceCpp(file.path(code.path, "EvaluationMetricsParallel.cpp"))

Cal_DEMaP <- function(seu, reduction = 'umap', ncores = min(10, RcppParallel::defaultNumThreads()-2), method = "C"){
  RcppParallel::setThreadOptions(numThreads = ncores)
  cat("calculate the graph and distance \n")
  edges <- as.Neighbor(seu@graphs$RNA_nn)@nn.idx
  edges <- as.vector(edges)
  edges <- data.frame("from" = rep(seq(1:ncol(seu)), 20),
                      "to" = edges,
                      "dist" = 1)
  graph <- cppRouting::makegraph(df = edges, directed = F)
  nodes <- unique(c(edges$from, edges$to))
  
  if (method == "R"){
    # geo_dist <- as.dist(cppRouting::get_distance_matrix(Graph = graph, from = nodes, to = nodes, allcores = T))
    # embed_dist <- dist(seu@reductions[[reduction]]@cell.embeddings)
    dist <- array(data = c(cppRouting::get_distance_matrix(Graph = graph, from = nodes, to = nodes, allcores = T), 
                           as.matrix(dist(seu@reductions[[reduction]]@cell.embeddings))), 
                  dim = c(ncol(seu), ncol(seu), 2))
    cat("calculate spearman correlation")
    spearman <- pbapply::pbapply(dist, 1, function(x){
      cor(x[, 1], x[, 2], method = "spearman")
    })
  }else if(method == "C"){
    from = nodes;to = nodes;
    from <- as.character(from)
    to <- as.character(to)
    allnodes <- c(from, to)
    from_id <- graph$dict$id[match(from, graph$dict$ref)]
    to_id <- graph$dict$id[match(to, graph$dict$ref)]
    if (ncores == 1){
      return(DEMaP(gfrom = graph$data[, 1], 
                   gto = graph$data[, 2], 
                   gw = graph$data[, 3],
                   NbNodes = graph$nbnod, 
                   dep = from_id, 
                   arr = to_id, 
                   embedding = seu@reductions[[reduction]]@cell.embeddings))
    }else{
      return(DEMaP_par(gfrom = graph$data[, 1], 
                       gto = graph$data[, 2], 
                       gw = graph$data[, 3],
                       NbNodes = graph$nbnod, 
                       dep = from_id, 
                       arr = to_id, 
                       embedding = seu@reductions[[reduction]]@cell.embeddings))
    }
  }
}

Cal_ARI <- function(seu, reduction = 'umap', method = "hierachy"){
  raw.cluster <- seu$seurat_clusters
  if (method == "hierachy"){
    embedding.cluster <- cutree(tree = hclust(dist(seu@reductions[[reduction]]@cell.embeddings)), 
                                k = length(levels(raw.cluster)))
  }else if(method == "kmeans"){
    embedding.cluster <- kmeans(x = seu@reductions[[reduction]]@cell.embeddings, 
                                centers = length(levels(raw.cluster)))$cluster
  }
  
  dict <- as.data.frame(table(raw.cluster, embedding.cluster))
  dict <- dplyr::arrange(dict, -dict$Freq)
  dict <- dict[!duplicated(dict$raw.cluster), ]
  aricode::ARI(c1 = raw.cluster, c2 = embedding.cluster)
}


