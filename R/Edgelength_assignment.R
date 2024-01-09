#' @title Assign barcode of unobserved node with parsimony
#'
#' @param phylo_c an object of class "phylo"
#'
#' @export
#'
#' @return a dataframe which stores the barcode of all nodes including observed leaves and unobserved internal node
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

#' @title Likelihood on cherry structure
#'
#' @param left_child n(length of site) x m(number of states in one site) binary matrix which saves the barcode of left node of a cherry
#' @param right_child n(length of site) x m(number of states in one site) binary matrix which saves the barcode of right node of a cherry
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param non_bifur_pro a parameter which describe the proportion of cell not bifurcated after one generation time.
#' @param branch_l branch length of left child
#' @param branch_r branch length of right child
#' @param r_height branch length incoming edges connecting the common ancestor node of two cherry to the fixed root
#'
#' @export
#'
#' @return a likelihood score on the cherry structure
#'
cherry_likelihood <- function(left_child,right_child,mu,alpha,non_bifur_pro,branch_l,branch_r,r_height){
  puni_constant <- log10(non_bifur_pro/(1-non_bifur_pro))
  parentnode_i <- get_parent(left_child,right_child,left_branch_length = branch_l,right_branch_length = branch_r,mu,alpha)
  site_num <- dim(left_child)[1]
  state_num <- dim(left_child)[2]
  parent_node_a <- array(data = 0,dim = c(site_num,state_num))
  for (i in 1:site_num)
    parent_node_a[i,] <- parentnode_i[i,] %*% t(get_probmatrix(t = r_height-min(branch_l,branch_r),mu = mu[i],alpha = alpha[[i]],state_num))
  score.a<- sum(log10(parent_node_a[,2])) + (branch_l+branch_r-2)*puni_constant
  return(score.a)
}

#' @title a bottom-up approach to assign the edge length under fixed tree height
#'
#' @param phylo an object of class "phylo"
#' @param nGen fixed tree height based on experimental duration and one generation of certain cell
#' @param state_num number of possible mutation outcomes
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param non_bifur_pro a parameter which describe the proportion of cell not bifurcated after one generation time.
#'
#' @export
#'
#' @description
#' a bottom-up approach to assign the edge length of a phylo object.The edge length will be the
#'     subtraction of vertex height of two node. At first, the vertex height of leaves are set
#'     as the fixed tree height, then the cherry structure is assigned an optimal
#'     tree height as the vertex height. We assume an peak in likelihood with fixed tree height.
#'
#' @return an object of class "phylo"

direct_assignment<- function(phylo,nGen,state_num,mu,alpha,non_bifur_pro){
  node_h <- list()
  n_sample <- length(phylo$tip.label)
  node_h$nheight[1:n_sample] <- nGen
  if(phylo$Nnode==n_sample) # 2023.01.08修正，未在github上更新
    phylo <- remove_pseudonode(phylo)
  node_path <- nodepath(phylo)
  node_df <- as.data.frame(sapply(node_path, "[", i = 1:max(sapply(node_path,length))))
  node_df[is.na(node_df)] <- 0
  barcode_state <- tidyr::unite(as.data.frame(get_node_character(phylo_c = phylo)),col = "state",sep = "")$state
  parent_node <- onehot_coding(prefix_state(node_info = barcode_state,state_num),state_num)
  puni_constant <- log10(non_bifur_pro/(1-non_bifur_pro))
  for(l in nrow(node_df):2)
    for (j in 1:(ncol(node_df)-1))
      for (k in (j+1):ncol(node_df))
        if (node_df[l,k] != 0 & (node_df[(l-1),j] == node_df[(l-1),k]) & (node_df[l,k] != node_df[l,j])){
          if (node_df[l,j] <= n_sample & node_df[l,k] <= n_sample){
            branch_l <- 0
          } else branch_l <- 1
          h_d <- node_h$nheight[node_df[l,j]]-node_h$nheight[node_df[l,k]] #  height difference in cell
          #limit <- min(node_h$nheight[c(node_df[l,j],node_df[l,k])])-1
          limit <- min(node_h$nheight[c(node_df[l,j],node_df[l,k])]) - (max(node.depth(phylo,method = 2))-node.depth(phylo,method = 2)[node_df[(l-1),k]]) # need further inspect
          while(branch_l < limit){
              score_c <- cherry_likelihood(left_child = parent_node[[node_df[l,j]]],right_child = parent_node[[node_df[l,k]]],branch_l = branch_l,branch_r = (branch_l+h_d),mu,alpha,non_bifur_pro,r_height = limit)
              score_n <- cherry_likelihood(left_child = parent_node[[node_df[l,j]]],right_child = parent_node[[node_df[l,k]]],branch_l = (1+branch_l),branch_r = (branch_l+h_d+1),mu,alpha,non_bifur_pro,r_height = limit)
              if (score_c > score_n){
                break
              } else branch_l <- branch_l + 1
            }
          node_h$nheight[node_df[(l-1),k]] <- min(node_h$nheight[c(node_df[l,j],node_df[l,k])]) - branch_l
      }
  for (x in 1:nrow(phylo$edge))
    phylo$edge.length[x] <- node_h$nheight[phylo$edge[x,2]]-node_h$nheight[phylo$edge[x,1]]
  phylo <- add_pseudonode(phylo)
  phylo$edge.length[1] <- node_h$nheight[n_sample+1]
  return(phylo)
}

#' Remove the pseudonode added
#'
#' @param phylo an object of "phylo"
#'
#' @return phylo
#' @export
remove_pseudonode <- function(phylo){
  n_sample <- length(phylo$tip.label)
  phylo$edge[which(phylo$edge[,1] > n_sample),1] <- phylo$edge[which(phylo$edge[,1] > n_sample),1] - 1
  phylo$edge[which(phylo$edge[,2] > n_sample),2] <- phylo$edge[which(phylo$edge[,2] > n_sample),2] - 1
  phylo$edge <- phylo$edge[-1,]
  if(length(phylo$edge.length) != 0)
    phylo$edge.length <- phylo$edge.length[-1]
  phylo$Nnode <-  phylo$Nnode-1
  return(phylo)
}
