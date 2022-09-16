library(igraph)
library(OmnipathR)
all_interactions <- import_all_interactions(
  organism = 9606,
  directed = 'no'
)

interactiongraph <- interaction_graph(all_interactions)
allshortestpaths <- all_shortest_paths(interactiongraph, c('WDR48'), c('PTEN'), mode = "ALL")$res
#res <-as.character(res)
allshortestpaths1 <- all_shortest_paths(interactiongraph, c('PTEN'), c('USP1'), mode = "ALL")$res
allshortestpaths1 <- allshortestpaths1[[1]]
allshortestpaths2 <- all_shortest_paths(interactiongraph, c('WDR48'), c('USP1'), mode = "ALL")$res
allshortestpaths2 <- c(allshortestpaths2[[1]], allshortestpaths2[[2]], allshortestpaths2[[3]])
allshortestpaths3 <- all_shortest_paths(interactiongraph, c('WDR48'), c('AKT1'), mode = "ALL")$res
allshortestpaths4 <- all_shortest_paths(interactiongraph, c('PTEN'), c('AKT1'), mode = "ALL")$res
allshortestpaths3 <- allshortestpaths3[[1]]
allshortestpaths4 <- allshortestpaths4[[1]]
allshortestpaths5 <- all_shortest_paths(interactiongraph, c('WDR48'), c('PHLPP1'), mode = "ALL")$res
allshortestpaths5 <- allshortestpaths5[[1]]
allshortestpaths7 <- all_shortest_paths(interactiongraph, c('PTEN'), c('PHLPP1'), mode = "ALL")$res
allshortestpaths7 <- c(allshortestpaths7[[21]], allshortestpaths7[[15]], allshortestpaths7[[11]])
res <- c(names(allshortestpaths), names(allshortestpaths1), names(allshortestpaths2), names(allshortestpaths3),
          names(allshortestpaths4), names(allshortestpaths5), names(allshortestpaths7))
edgelist2 <- matrix(c(res[4], res[5], res[5], res[1], res[2], res[3], res[2], res[26], 
                      res[1], res[15], res[4], res[14], res[14], res[15], res[27], res[14], res[14], res[29], res[15], res[29]), nc = 2, byrow = TRUE)
graph2 <-graph_from_edgelist(edgelist2)
plot(as.undirected(graph2))

