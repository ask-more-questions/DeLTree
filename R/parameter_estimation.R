#' Assign barcode with Parsimony
#'
#' @param phylo_c a object of "phylo"
#'
#' @export
#'
#' @return a dataframe which stores the the barcode of leaves and nodes
get_node_character <- function(phylo_c){
  tlabels <- sapply(sapply(phylo_c$tip.label,function(x){strsplit(x,split = "_")}),"[")[2,]
  tip_character <- t(sapply(strsplit(tlabels,split = ""),"["))
  npath <- nodepath(phylo_c)
  node_order <- seq((length(phylo_c$tip.label)+1),length.out = phylo_c$Nnode)
  node_character <- c()
  for( i in 2:length(node_order))
  {node_to_tpindex <- which(sapply(sapply(npath,function(x){x %in% node_order[i]},simplify = FALSE),function(x){TRUE %in% x}))
  tip_info <- sapply(strsplit(tlabels[node_to_tpindex],split = ""),"[")
  n_mtx <- apply(tip_info,1,function(x){
    if(length(unique(x))==1) return(unique(x)) else return("1")})
  node_character <- rbind(node_character,n_mtx)
  }
  # add root
  node_character <- rbind(rep(1,ncol(tip_character)),node_character)
  tree_character <- rbind(tip_character,node_character)
  rownames(tree_character) <- 1:max(node_order)
  return(tree_character)
}


#' Count the number of mutational events
#'
#' @param phylo_c an object of class "phylo"
#' @param site_num the length of barcode
#' @param state_num  number of possible mutation outcomes
#'
#' @export
#'
#' @return a matrix which stores the number of mutational/non mutation events on edge of the tree

get_site_mutation <- function (phylo_c, site_num, state_num) {
  tree_character <- get_node_character(phylo_c)
  s <- apply(tree_character, 2, as.numeric)
  n_sample <- length(phylo_c$tip.label)
  one_generation_edge <- which(phylo_c$edge.length==1)
  count_info <- matrix(data = 0, nrow = state_num, ncol = site_num)
  if(length(one_generation_edge)!=0){
     parent_child_mutation <- matrix(data = 0, nrow = length(one_generation_edge),ncol = site_num)
  for (site in 1:site_num) {
    for (i in 1:length(one_generation_edge)) {
      node1 <- phylo_c$edge[one_generation_edge[i], 1]
      node2 <- phylo_c$edge[one_generation_edge[i], 2]
      parent_child_mutation[i, site] <- s[node1, site] - s[node2, site]
    }
    for (j in 1:state_num) {
      count_info[j, site] <- length(which(parent_child_mutation[,site] == (1 - j + 1)))
    }
    count_info[2,site] <- sum(s[phylo_c$edge[one_generation_edge,2], site]==1)#only consider edge length 1 situation
  }
  }
  return(count_info)
}
