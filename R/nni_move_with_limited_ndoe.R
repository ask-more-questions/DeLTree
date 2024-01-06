# nni move for tree_topology
#' @title Identify node within the replicate cluster
#'
#' @param phylo  an object of class "phylo"
#'
#' @export
#'
#' @return a vector which stores the node within a subtree full of replicate barcode sequence
within_replicate_node <- function(phylo){
  n_sample <- length(phylo$tip.label)
  internal_node <- (n_sample+2):max(phylo$edge[,1])
  within_replicate_node <- c()
  for (i in internal_node){
    child_leaves <- Descendants(phylo,node = i,type = "tips")[[1]]
    child_code <- unique(sapply(strsplit(phylo$tip.label[child_leaves],split = "_"),"[",2))
    if (length(child_code) == 1)
      within_replicate_node <- append(within_replicate_node,values = i)
  }
  return(within_replicate_node)
}




#' Selet nni node
#'
#' @param phylo an object of class "phylo"
#'
#' @return a list of nni trees
#' @export
get_nni_node_noRT <- function(phylo){
  edgelength_c <- cbind(phylo$edge,phylo$edge.length)
  within_replicate_node <- within_replicate_node(phylo)
  nni_edge <- edgelength_c[edgelength_c[,2] > length(phylo$tip.label),]
  nni_edge <- nni_edge[!(nni_edge[,1] %in% within_replicate_node),]
  nni_edge <- nni_edge[-1,]
  nni_node <- matrix(data = 0, nrow = 2*nrow(nni_edge),ncol = 5)
  j = 1
  for (i in 1:nrow(nni_edge)){
    p_node <- nni_edge[i,1];nni_node[c(j,j+1),1] <- p_node #
    c1_node <- nni_edge[i,2];nni_node[c(j,j+1),2] <- c1_node
    c2_node  <- setdiff(edgelength_c[(edgelength_c[,1] %in% p_node),2],c1_node);nni_node[c(j,j+1),3] <- c2_node #
    cc_node_index <- which(edgelength_c[,1] %in% c1_node)
    c11_node <- edgelength_c[cc_node_index[1],2];nni_node[j,4] <- c11_node;nni_node[j+1,5] <- c11_node
    c12_node <- edgelength_c[cc_node_index[2],2];nni_node[j,5] <- c12_node;nni_node[j+1,4] <- c12_node
    j<-j+2
  }
  return(nni_node)
}

#' @title Tree rearrangements of one NNI(nearest neighbor interchange) move without node within certain cluster
#'
#' @param phylo  an object of class "phylo".
#'
#' @export
#'
#' @return a list of nni trees
nni_tree_noRT <- function(phylo){
  nni_node <- get_nni_node_noRT(phylo)
  nni_treelist <- list()
  for (i in 1:nrow(nni_node)){
    nni_treelist[[i]] <- phylo
    edgelength_c <- cbind(phylo$edge,phylo$edge.length)
    p_node <- nni_node[i,1]
    c1_node <- nni_node[i,2]
    l <- edgelength_c[edgelength_c[,2] %in% c1_node,3]
    edgelength_c[edgelength_c[,2] %in% c1_node,3] <- 0
    c2_node <- nni_node[i,3]
    c11_node <- nni_node[i,4]
    edgelength_c[edgelength_c[,2] %in% c2_node,1] <- c1_node # edgelength will not change
    edgelength_c[edgelength_c[,2] %in% c11_node,1] <- p_node # here this node of the original cherry is outgrouped
    c11_end_node <- unlist(sapply(nodepath(phylo),function(x){if(c11_node %in% x) tail(x,1)}))
    edgelength_c[(edgelength_c[,2] %in% c11_end_node),3] <- edgelength_c[(edgelength_c[,2] %in% c11_end_node),3] + l
    c12_node <- nni_node[i,5]
    edgelength_c[edgelength_c[,2] %in% c12_node,1] <- c1_node
    c12_end_node <- unlist(sapply(nodepath(phylo),function(x){if(c12_node %in% x) tail(x,1)}))
    edgelength_c[(edgelength_c[,2] %in% c12_end_node),3] <- edgelength_c[(edgelength_c[,2] %in% c12_end_node),3] + l
    nni_treelist[[i]]$edge <- edgelength_c[,c(1,2)]
    nni_treelist[[i]]$edge.length <- edgelength_c[,3]
  }
  return(nni_treelist)
}


#' @title Likelihood score
#'
#' @param phylo  an object of class "phylo".
#' @param parent.node a list of n(length of site) x m(number of states in one site) binary matrix which saves the barcode of leaves
#' @param node_edge_length a data frame which bind the edge attribute and edge.length attribute of a phylo structure.
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#'
#' @export
#'
#' @return A list of matrix of number of nodes

emurate_node_path_withpseudonode <- function(phylo,parent.node,node_edge_length,mu,alpha){
  node_path <- nodepath(phylo)
  node_df <- as.data.frame(sapply(node_path, "[", i = 1:max(sapply(node_path,length))))
  site_num <- dim(parent_node[[1]])[1]
  state_num <- dim(parent_node[[1]])[2]
  node_df[is.na(node_df[,])] <- 0# drop na value as 0
  node_dep <- node.depth(phylo,method = 2)
  node_max <- max(node_dep)
  n_sample <- length(node_path)
  ngen <- (nrow(node_df)-1)
  # core function
  for (g in 1:ngen){
    current_g <- (ngen - g + 2)
    for (i in 1:(n_sample-1))
      for (j in (i+1):n_sample){
        criteria_1 <- (node_df[(current_g-1),i] == node_df[(current_g-1),j])
        criteria_2 <- (node_df[(current_g),i] == node_df[(current_g),j])
        criteria_3 <- (node_df[(current_g),i] > 0) # mark here : change >= 0 to >0
        if (criteria_1 & !criteria_2 & criteria_3)
        {left_index <- node_df[(current_g),i]
        right_index <- node_df[(current_g),j]
        parent_index  <- node_df[(current_g-1),i]
        br_left_index <- (node_edge_length[,1]  %in%  parent_index) & (node_edge_length[,2] %in% left_index)
        br_right_index <- (node_edge_length[,1]  %in% parent_index) & (node_edge_length[,2] %in% right_index)
        parent_node[[parent_index]] <- get_parent(parent_node[[left_index]],parent_node[[right_index]],
                                                  node_edge_length[br_left_index,3], node_edge_length[br_right_index,3],mu,alpha)
        }
      }
  }
  parent_node[[n_sample+1]] <- matrix(data = 0,nrow = site_num,ncol = state_num)
  for (i in 1:site_num)
    parent_node[[n_sample+1]][i,] <- parent_node[[n_sample+2]][i,] %*% t(get_probmatrix(t = node_edge_length[1,3],mu = mu[i],alpha = alpha[[i]],state_num))
  return(parent_node)
}

#' @title Initial edgelength assignment with longest pending edge
#'
#' @param phylo an object of class "phylo".
#' @param nGen fixed tree height based on experimental duration and one generation of certain cell
#'
#' @export
#'
#' @return a data frame which bind the edge attribute and edge.length attribute of a phylo structure.

initial_edgelength_pseudonode <- function(phylo,nGen){
  nsample <-length(phylo$tip.label)
  phylo$edge.length <- rep(1,nrow(phylo$edge))
  taxa_index <- which(phylo$edge[,2] %in% 1:nsample)
  pending_edgelength <- sapply(nodepath(phylo),function(x){return(nGen+3-length(x))})
  phylo$edge.length[taxa_index] <- pending_edgelength
  edgelength <- cbind(phylo$edge,phylo$edge.length)
  edgelength[1,3] <- 0
  return(edgelength)
}

#' @title Find non-bifurcation node
#'
#' @param edgelength a data frame which bind the edge attribute and edge.length attribute of a phylo structure.
#'
#' @export
#'
#' @return a vector of node index

get_sparenode_name <- function(edgelength){
  free_edge_index <- which(edgelength[,3] > 1 | (edgelength[,2] <= edgelength[1,1] & edgelength[,3] == 1))
  free_node_name <- names(which(table(edgelength[free_edge_index,1])==2))
  return(free_node_name)
}

#' @title Bottom up edge length iteration which reduce non bifurcation event
#'
#' @param edgelength a data frame which bind the edge attribute and edge.length attribute of a phylo structure.
#' @param free_node_initial  index of parent node whose two descendant node both have non bifurcation events
#'
#' @export
#'
#' @return a data frame which bind the edge attribute and edge.length attribute of a phylo structure.
edgelength_inter <- function(edgelength,free_node_initial){
  edgelength_c <- edgelength
  indegree_1p_index <- which(edgelength_c[,2] %in% free_node_initial) # len 1
  outdegree_1m_index <- which(edgelength_c[,1] %in% free_node_initial) # len 2
  pending_edge_1 <- (sum(edgelength_c[outdegree_1m_index,2] < edgelength_c[1,1])==2 & sum(edgelength_c[outdegree_1m_index,3] ==1)==2 )
  if (sum(edgelength_c[outdegree_1m_index,3] > 1)==2 | pending_edge_1){
    edgelength_c[indegree_1p_index,3] <- edgelength_c[indegree_1p_index,3] + 1
    edgelength_c[outdegree_1m_index,3] <- edgelength_c[outdegree_1m_index,3] -1
  }
  return(edgelength_c)
}
#' @title Bottom up edgelength optimzation
#'
#' @param tree an object of class "phylo"
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param nGen fixed tree height based on experimental duration and one generation of certain cell
#' @param non_bifur_pro A parameter which describes the proportion of cell not bifurcated after one generation time
#' @param state_num number of mutation outcomes
#'
#' @export
#'
#' @return a data frame which bind the edge attribute and edge.length attribute of a phylo structure.

local_optiaml_withpseudoRT <- function(tree,mu,alpha,nGen,non_bifur_pro,state_num){
  state <- sapply(strsplit(x=tree$tip.label,split = "_"),"[")[2,]
  node_info <- t(sapply(strsplit(state,split = ""),"["))
  parent.node <- onehot_coding(prefix_state(node_info,state_num),state_num)
  n_sample <- length(tree$tip.label)
  edgelength <- initial_edgelength_pseudonode(tree,nGen)
  tree.freenode <- sort(get_sparenode_name(edgelength),decreasing = TRUE)
  free_node_num <- length(tree.freenode)
  edgelength_c <-  edgelength
  initial_node <- tree.freenode
  move <- 0
  result <- list()
  while(length(initial_node!=0) & initial_node[1] != (n_sample+1)){
    move <- move + 1
    likelihood_c <- emurate_node_path_withpseudonode(phylo = tree,parent.node,node_edge_length = edgelength_c,mu = mu,alpha = alpha)
    score_c <- sum(log10(likelihood_c[[n_sample+1]][,2])) + log10(bifur_punish(bifurcation_pro = non_bifur_pro,edgelength = edgelength_c))
    new_edgelength <- edgelength_inter(edgelength = edgelength_c,free_node_initial = initial_node[1])
    new_likelihood <- emurate_node_path_withpseudonode(phylo = tree,parent.node,node_edge_length = new_edgelength,mu = mu,alpha = alpha)
    score_new <- sum(log10(new_likelihood[[n_sample+1]][,2])) + log10(bifur_punish(bifurcation_pro = non_bifur_pro,edgelength = new_edgelength))
    c <- 1
    while ((c <= nGen)&(score_new > score_c))
    {edgelength_c <- new_edgelength
    c <- c+1
    likelihood_c <- emurate_node_path_withpseudonode(phylo = tree,parent.node,node_edge_length  = edgelength_c,mu = mu,alpha = alpha)
    score_c <- sum(log10(likelihood_c[[n_sample+1]][,2])) + log10(bifur_punish(bifurcation_pro = non_bifur_pro,edgelength = edgelength_c))
    new_edgelength <- edgelength_inter(edgelength = edgelength_c,free_node_initial = initial_node[1])
    new_likelihood <- emurate_node_path_withpseudonode(phylo = tree,parent.node,node_edge_length = new_edgelength,mu = mu,alpha = alpha)
    score_new <- sum(log10(new_likelihood[[n_sample+1]][,2])) + log10(bifur_punish(bifurcation_pro = non_bifur_pro,edgelength = new_edgelength))
    }
    tree.freenode
    new_node <- sort(get_sparenode_name(edgelength_c),decreasing = TRUE)
    new_node <- setdiff(new_node,tree.freenode)
    tree.freenode <- append(tree.freenode,new_node)
    if(length(new_node) != 0)
      initial_node[1] <- new_node
    else
      initial_node <- setdiff(initial_node,initial_node[1])
  }
  return(edgelength_c)
}

#' @title NNI move and edge length optimization
#'
#' @param current_tree an object of class "phylo".
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param nGen fixed tree height based on experimental duration and one generation of certain cell
#' @param non_bifur_pro A parameter which describes the proportion of cell not bifurcated after one generation time
#' @param state_num number of mutation outcomes
#'
#' @export
#'
#' @return a list of nni trees with optimized edge length and likelihood.

nni_iter_withedgelength_pseudonode <- function(current_tree,mu,alpha,nGen,non_bifur_pro,state_num){
  nni_recorder <- list()
  state <- sapply(strsplit(x=current_tree$tip.label,split = "_"),"[")[2,]
  node_info <- t(sapply(strsplit(state,split = ""),"["))
  parent.node <- onehot_coding(prefix_state(node_info,state_num),state_num)
  n_sample <- length(parent.node)
  nni_tree <-  nni_tree_noRT(current_tree)
  # node_depth <- sapply(nni_tree,function(x){max(node.depth(x,method = 2))-1})
  # nni_tree <- nni_tree[node_depth <= nGen]
  nni_likelihood <- c()
  for (i in 1:length(nni_tree)){
    edgelength_c <- local_optiaml_withpseudoRT(tree=nni_tree[[i]],mu = mu,alpha=alpha,nGen,non_bifur_pro,state_num)
    nni_tree[[i]]$edge <- edgelength_c[,1:2]
    nni_tree[[i]]$edge.length <- edgelength_c[,3]
    nni_likelihood[i] <- sum(log10(emurate_node_path_withpseudonode(phylo = nni_tree[[i]],parent.node,node_edge_length = edgelength_c,mu,alpha)[[n_sample+1]][,2])) + log10(bifur_punish(bifurcation_pro = non_bifur_pro,edgelength = edgelength_c))
  }
  nni_best <- which(nni_likelihood == max(nni_likelihood))[1]
  nni_best_tree <- nni_tree[[nni_best]]
  nni_recorder$likelihood <- nni_likelihood
  nni_recorder$best_tree <- nni_best_tree
  return(nni_recorder)
}
