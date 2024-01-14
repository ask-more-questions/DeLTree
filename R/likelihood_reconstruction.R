#' Irreversible consensus
#'
#' @param node_info  A n(cell_num) x m(site_num) matrix which stores the barcode
#' @param site_index  The index of two descendant cell of a node
#' @param site_num The length of barcode
#' @description
#' This function applies a parsimony method to assign the barcode with the irreversible rule of mutation and the barcode of its two direct descendant.
#'
#' @export
#'
#' @return A vector of with the length of site_num which stores the barcode of the parent node
get_consensus <- function(node_info,site_index,site_num){
  parent <- c()
  node1 <- node_info[site_index[1],]
  node2 <- node_info[site_index[2],]
  site <- abs(node1-node2)
  for (s in 1:site_num){
    if (site[s] == 0) parent[s] <- node1[s] else parent[s] <- 1
  }
  return(parent)
}


#' One hot encoding
#'
#' @param processed_tip_label barcode sequence with a prefix
#' @param state_num number of mutation outcomes
#' @description
#'     By default, this function take a matrix which stores the barcode information of observed cell and convert the barcode sequence into m(length of the barcode or mutation site) x n(number of mutation outcomes) binary matrix with one hot encoding techniques,
#'     the binary matrix of barcode are saved in the list of n cells.
#' @export
#'
#' @return A list with the length of cell number. each element include a site_num x state_num binary matrix

onehot_coding<- function(processed_tip_label,state_num){
  tip.label_list <- strsplit(processed_tip_label,split = "")
  n_sample <- length(processed_tip_label)
  parent.node <-list()
  for (i in 1:n_sample){
    dummy <- caret::dummyVars("~.",data = tip.label_list[i])
    parent.node[[i]] <- tail(predict(dummy,newdata = tip.label_list[i]),-state_num)

  }
  return(parent.node)
}
#' Title
#'
#' @param node_info A n(cell_num) x m(site_num) matrix which stores the barcode
#' @param state_num number of mutation outcomes
#' @description
#' A preprocess step which adds a prefix which contains all possible mutation outcomes to ensure the dimention of binary matrix before one hot encoding
#'
#' @export
#'
#' @return the barcode sequence with a prefix which contains all possible mutation outcomes
prefix_state <- function(node_info,state_num){
  state <- tidyr::unite(as.data.frame(node_info),col = "state",sep = "")
  prefix <- paste0(c(0:(state_num-1)),collapse = "")
  state_seq <- paste(prefix,state$state,sep = "")
  return(state_seq)
}



#' @title Probability matrix of dicrete edge length
#'
#' @param t integer.number of generation or number of bifurcation
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param state_num number of mutation outcomes.
#'
#' @export
#'
#' @return A probability matrix

get_probmatrix <- function(t,mu,alpha,state_num){
  if(length(alpha) != 1){
    if(round(sum(alpha),digits = 3) !=1)
      stop("In multiple sites : vector alpha should sum up to 1")
    p_matrix <- diag(1,nrow = state_num,ncol = state_num)
    p_matrix[2,] <- append(x = mu*alpha,values = (1-mu),after = 1)
  } else p_matrix <- matrix(c(1,0,0,mu*alpha,(1-mu),mu*(1-alpha),0,0,1),nrow = 3,byrow = TRUE)
  if(round(sum(p_matrix[2,]),digits = 3)!=1){
    stop("Wrong P_matrix : rowsum not 1")
  }
  egi_value <- eigen(p_matrix)$values
  egi_vector <- eigen(p_matrix)$vectors
  egi_vector_inver <- solve(egi_vector,diag(x=1,ncol = state_num,nrow = state_num))
  p_n <- egi_vector %*% diag(egi_value^t,nrow = state_num,ncol = state_num) %*% egi_vector_inver
  return(p_n)
}


#' Likelihood computation with matrix
#'
#' @param child_left n(length of site) x m(number of states in one site) binary matrix which saves the barcode of left node of a cherry
#' @param child_right n(length of site) x m(number of states in one site) binary matrix which saves the barcode of right node of a cherry
#' @param left_branch_length left branch length
#' @param right_branch_length right branch length
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#'
#' @return n(length of site) x m(number of states in one site) binary matrix which saves the barcode of common ancestor node of a cherry
#' @export
get_parent <- function(child_left,child_right,left_branch_length,right_branch_length,mu,alpha){
  site_num <- dim(child_left)[1]
  state_num <- dim(child_left)[2]
  parent <- array(data = 0,dim = c(site_num,state_num))
  constant.pr <- diag(1,nrow = state_num,ncol = state_num)
  for (i in 1:site_num){
    if (left_branch_length > 0) left_tree <- child_left[i,] %*% t(get_probmatrix(t = left_branch_length,mu = mu[i],alpha = alpha[[i]],state_num)) else left_tree <- child_left[i,] %*% constant.pr
    if (right_branch_length > 0) right_tree <- child_right[i,] %*%  t(get_probmatrix(t= right_branch_length,mu = mu[i], alpha = alpha[[i]],state_num)) else right_tree <- child_right[i,] %*% constant.pr
    parent[i,] <- left_tree * right_tree
  }
  return(parent)
}

#' @title Topology biased score
#'
#' @param non_bifur_pro a parameter which describes the proportion of cell not bifurcated after one generation time.
#' @param edgelength a data frame which bind the edge attribute and edge.length attribute of a phylo structure.
#' @description
#'    This function takes the proportion of  bifurcation event/non bifurcation event
#'    as the base and number of non-bifurcation event as the power to propose a topology biased
#'    score under the discrete edge length and fixed tree height model in which
#'    edge length bigger than one implies (n-1) continuous non bifurcation event.
#' @export
#'
#' @return a length-one numeric
bifur_punish<- function(non_bifur_pro,edgelength){
  edgelength <- edgelength[-1,]
  non_bifurcation_time <- sum(edgelength[(edgelength[,3]>1),3]-1)
  #reward <- sum(edgelength[,3] == 0)
  y <- (non_bifur_pro/(1-non_bifur_pro))^non_bifurcation_time
  return(y)
}

#' @title Likelihood of bifurcation within two generations
#'
#' @param left_child n(length of site) x m(number of states in one site) binary matrix which saves the barcode of left node of a cherry
#' @param right_child n(length of site) x m(number of states in one site) binary matrix which saves the barcode of right node of a cherry
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param non_bifur_pro a parameter which describe the proportion of cell did not bifurcation after one generation time.
#'
#' @description
#' this function compute a pairwise likelihood of two cells orgins from a common ancestor within one bifurcation time. Given a tree height 2,
#'     there is an optional unobserved ancestor assiociated with an bifurcation event in the cell. The subtraction of having and
#'     not having an internal node are used as the pairwise distance of two barocdes in two generation time.
#'
#' @import ape
#' @import phangorn
#' @importFrom stats predict runif
#' @importFrom utils tail
#'
#' @export
#'
#' @return a non-negative score which describe  how many times its more likely to bifurcate in one generation time than not to bifurcate.
str_cherry_lik <- function(left_child,right_child,mu,alpha,non_bifur_pro){
  puni_constant <- log10(non_bifur_pro/(1-non_bifur_pro))
  site_num <- dim(left_child)[1]
  state_num <- dim(left_child)[2]
  parentnode_i <- get_parent(left_child,right_child,left_branch_length = 2,right_branch_length = 2,mu,alpha)
  parent_node_c <- array(data = 0,dim = c(site_num,state_num))
  for (i in 1:site_num)
    parent_node_c[i,] <- parentnode_i[i,] %*% t(get_probmatrix(t = 0,mu = mu[i],alpha = alpha[[i]],state_num))
  score.c <- sum(log10(parent_node_c[,2]))+puni_constant
  parentnode_i <- get_parent(left_child,right_child,left_branch_length = 1,right_branch_length = 1,mu,alpha)
  parent_node_a <- array(data = 0,dim = c(site_num,state_num))
  for (i in 1:site_num)
    parent_node_a[i,] <- parentnode_i[i,] %*% t(get_probmatrix(t = 1,mu = mu[i],alpha = alpha[[i]],state_num))
  score.a <- sum(log10(parent_node_a[,2]))
  parentnode_i <- get_parent(left_child,right_child,left_branch_length = 0,right_branch_length = 0,mu,alpha)
  parent_node_b <- array(data = 0,dim = c(site_num,state_num))
  for (i in 1:site_num)
    parent_node_b[i,] <- parentnode_i[i,] %*% t(get_probmatrix(t = 2,mu = mu[i],alpha = alpha[[i]],state_num))
  score.b <- sum(log10(parent_node_b[,2]))
  score_final <- max(score.a,score.b,score.c)-score.c
  return(score_final)
}


#' Add pseudo node
#'
#' @param phylo an object of class "phylo".
#'
#' @description
#' A function which adds a node to the root which allows flexible observation time
#'
#' @export
#'
#' @return a phylo structure


add_pseudonode <- function(phylo){
  n_sample <- length(phylo$tip.label)
  phylo$edge[which(phylo$edge[,1] > n_sample),1] <- phylo$edge[which(phylo$edge[,1] > n_sample),1] + 1
  phylo$edge[which(phylo$edge[,2] > n_sample),2] <- phylo$edge[which(phylo$edge[,2] > n_sample),2] + 1
  phylo$edge <- rbind(c(n_sample+1,n_sample+2),phylo$edge)
  if(length(phylo$edge.length) != 0)
    phylo$edge.length <- append(phylo$edge.length,values = 0,after = 0)
  phylo$Nnode <-  n_sample
  return(phylo)
}




#' Tree Reconstruction
#'
#' @param character_info A data frame with two columns which saves cell name and barcode respectively
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param non_bifur_pro  A parameter which describes the proportion of cell not bifurcated after one generation time
#' @param alternative_threshold Default as zero, if 1, two sequences are 10 times more likely to have a bifurcation event in 2 generation.
#' @param nGen fixed tree height based on experimental duration and one generation of certain cell
#' @param site_num  the length of barcode
#' @param state_num number of mutation outcomes
#' @param precluster_replicate Logical. whether or not to cluster the cells with identical barcodes as cluster before reconstruct phylogeny.default as TRUE.
#'
#' @description
#' This function reconstruct the tree based on the pairwise likelihood of having one bifurcation in two generation time under the bottom up approach.
#'     The barcode of unobserved internal vertex follows the irreversible mutation rule and the fixed root state.
#' @import caret
#' @import phytools
#' @export
#'
#' @return a phylo structure without edge length

likelihood_based_recon <- function(character_info,mu,alpha,non_bifur_pro,alternative_threshold,nGen,site_num,state_num,precluster_replicate = TRUE){
  node_info <- apply(sapply(strsplit(character_info$state,split = ""),"["),1,as.integer)
  n_node <- nrow(character_info)
  recorder <- list()
  if (precluster_replicate){
    for(i in 1:n_node){
      recorder$cherry_index[[i]] <- paste(character_info$cell[i],character_info$state[i],sep = "_")
      recorder$cherry_code[[i]] <- node_info[i,]
      recorder$round[[i]] <- 0
    }
    unique_state <- names(table(character_info$state))
    for (i in 1:length(unique_state)){
      recorder$cherry_code[[n_node+i]] <- as.integer(strsplit(unique_state[i],split = "")[[1]])
      replicate_group <- paste((character_info$cell[which(character_info$state %in% unique_state[i])]),unique_state[i],sep = "_")
      recorder$round[[n_node+i]] <- ceiling(log2(length(replicate_group)) + 1)
      while(length(replicate_group)>=4){
        for (t in 1:floor(length(replicate_group)/2)){
          replicate_group[t] <- paste("(",paste(replicate_group[t],replicate_group[t+1],sep = ","),")",sep = "")
          replicate_group <- replicate_group[-(t+1)]
        }
      }
      if (length(replicate_group)==3){
        replicate_group <- paste("(",paste("(",paste(replicate_group[1],replicate_group[2],sep = ","),")",sep = ""),",",replicate_group[3],")",sep = "")
      }
      if(length(replicate_group)==2){
        replicate_group <- paste("(",paste(replicate_group[1],replicate_group[2],sep = ","),")",sep = "")
      }
      recorder$cherry_index[[n_node+i]] <- replicate_group
    }
    ntaxa <- length(unique_state)
    n_node <- n_node+ntaxa
    round <- max(1,min(unlist(sapply(recorder$round,"["))[1:n_node][(unlist(sapply(recorder$round,"["))[1:n_node]) >0]))
    # start the round
  } else {
    for(i in 1:n_node){
      recorder$cherry_index[[i]] <- paste(character_info$cell[i],character_info$state[i],sep = "_")
      recorder$cherry_code[[i]] <- node_info[i,]
      recorder$round[[i]] <- 1
    }
    ntaxa <- n_node
    round <- 1
  }
  round_limit <- max(unlist(recorder$round))
  while(ntaxa != 1 | round <= round_limit){
    round_index <- which(unlist(sapply(recorder$round,"["))[1:n_node] == round)
    while(length(round_index) <= 1){
      recorder$round[[round_index]] <- recorder$round[[round_index]]+1
      round <- round + 1
      round_index <- which(unlist(sapply(recorder$round,"["))[1:n_node] == round)
    }
    node_info <- t(sapply(recorder$cherry_code,"["))[round_index,]
    round <- round + 1
    ntaxa <- nrow(node_info)
    # mark here
    dist_pairwise <- matrix(data=0, ncol = ntaxa,nrow = ntaxa)
    parent_node <- onehot_coding(processed_tip_label = prefix_state(node_info,state_num),state_num)
    for (i in 1:(ntaxa-1)){
      for (j in (i+1):ntaxa){
        left_child <- parent_node[[i]]
        right_child <- parent_node[[j]]
        dist_pairwise[i,j] <- str_cherry_lik(left_child,right_child,mu,alpha,non_bifur_pro)
      }
    }
    left_taxa <- c(1:ntaxa)
    diag(dist_pairwise) <- -1
    n=1
    new_threshold <- min(min(dist_pairwise[dist_pairwise > 0]),alternative_threshold)
    while (max(dist_pairwise) > new_threshold | (n==1)){
      ttt1 <- as.vector(which(dist_pairwise == max(dist_pairwise),arr.ind = TRUE)[1,])
      recorder$cherry_index[[n_node+n]] <- paste("(",paste(recorder$cherry_index[[round_index[ttt1[1]]]],recorder$cherry_index[[round_index[ttt1[2]]]],sep = ","),")",sep = "")
      recorder$cherry_code[[n_node+n]] <- get_consensus(node_info = node_info,site_index = ttt1,site_num)
      recorder$round[[n_node+n]] <- round
      dist_pairwise[ttt1,] <- 0
      dist_pairwise[,ttt1] <- 0
      left_taxa <- setdiff(left_taxa,ttt1)
      n=n+1
    }
    if (length(left_taxa) != 0){
      for (i in 1:length(left_taxa)){
        recorder$cherry_index[[n_node+n]] <- recorder$cherry_index[[round_index[left_taxa[i]]]] #add n_node-ntaxa
        recorder$cherry_code[[n_node+n]] <- recorder$cherry_code[[round_index[left_taxa[i]]]]# add n_node-ntaxa
        recorder$round[[n_node+n]] <- round
        n=n+1
      }
    }
    n=n-1
    # node_info <- t(sapply(recorder$cherry_code,"["))[(n_node+1):(n_node+n),]
    n_node <- n_node+n
    ntaxa <- length(which(unlist(sapply(recorder$round,"["))[1:n_node] == round))
  }
  tree_txt <- paste(recorder$cherry_index[[n_node]],";root",sep = "")
  recon_tree <- read.newick(text = tree_txt)
  return(recon_tree)
}




