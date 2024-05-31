

#' @title a forward approach to simulate tree topology with bifurcation
#'
#' @param bifur_pro a paramter which describe the proportion of cell bifurcated after one generation time.
#' @param nGen fixed tree height of n generation time
#'
#' @return an object of class "phylo"
#'
#' @export
topology_simulation <- function(bifur_pro,nGen){
  tree_i <- stree(n = 1,tip.label = 1)
  tree_i$edge.length <- rep(0,nrow(tree_i$edge))
  round <- 0
  add_n <- 1
  while(round < nGen){
    leaves <- tree_i$tip.label
    for(i in 1:length(leaves)){
      l <- runif(1,min = 0,max = 1)
      if (l <= bifur_pro){
        tree_i <- add.tips(tree_i,tips = c((add_n+1),(add_n+2)),where = as.character(leaves[i]),edge.length = c(1,1))
        tree_i <- drop.tip(tree_i,tip = as.character(leaves[i]),collapse.singles = FALSE)
        add_n <- add_n + 2
      } else {
        leaf_index <- which(tree_i$tip.label %in% as.character(leaves[i]))
        tree_i$edge.length[which(tree_i$edge[,2] %in% leaf_index)] <- tree_i$edge.length[which(tree_i$edge[,2] %in% leaf_index)] + 1
      }
    }
    round <- round+1
  }
  return(tree_i)
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



#' @title Lineage barcode simulation
#'
#' @param branch_length branch length
#' @param parent_state n(length of site) x m(number of states in one site) binary matrix which saves the barcode of parent node
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#' @param site_num the length of barcode
#' @param state_num number of possible mutation outcomes
#'
#' @return a binary matrix which stores the barcode
#'
#' @export
node_state_sim <- function(branch_length,parent_state,mu,alpha,site_num,state_num){
  node_state <- matrix(data = rep(rep(0,state_num),each = site_num),nrow =site_num,byrow = FALSE)
  yesorno <- runif(site_num,0,1)
  #intermediate_state <-  matrix(data = rep(rep(0,state_num),each = site_num),nrow =site_num,byrow = FALSE)
  for (i in 1:site_num){
    intermediate_state <- parent_state[i,] %*% get_probmatrix(t = branch_length,mu[i],alpha[[i]],state_num)
    condition <- sapply(1:state_num,function(x){return(sum(intermediate_state[1:x]))})
    s_index <-  which(condition-yesorno[i] >0)[1]
    node_state[i,s_index] <- 1
  }
  return(node_state)
}


#' @title reverse onehot encoding
#'
#' @param n_state_matrix n(length of site) x m(number of states in one site) binary matrix which saves the barcode
#' @param state_num number of possible mutation outcomes
#'
#' @return a vector of barcode
#'
#' @export

reverse_onehot_encoding  <-  function(n_state_matrix,state_num){
  state <- c()
  for (i in 1:nrow(n_state_matrix)){
    s_index <- which(n_state_matrix[i,] == 1)
    state[i] <- s_index-1
  }
  return(state)
}

#' Lineage simulation
#'
#' @param tree an object of class "phylo".
#' @param state_num number of possible mutation outcomes
#' @param site_num the length of barcode
#' @param mu a vector of site specific mutation probability.
#' @param alpha a list of vectors which describe the site specific priors of mutation outcomes.
#'
#' @return a data frame with two columns which saves cell name and barcode of leaves respectively
#'
#' @export

lineage_sim <- function(tree,state_num,site_num,mu,alpha){
  n_cell <- length(tree$tip.label)
  node_sim_list <- list()
  node_edgelength <- cbind(tree$edge,tree$edge.length)
  rt_index <- node_edgelength[1,1]
  node_sim_list[[rt_index]] <- matrix(data = rep(append(rep(0,(state_num-1)),values = 1,after = 1),each = site_num),nrow =site_num,byrow = FALSE)
  for (i in 1:nrow(node_edgelength)){
    parent_index <- node_edgelength[i,1]
    child_index <- node_edgelength[i,2]
    node_sim_list[[child_index]] <- node_state_sim(branch_length = node_edgelength[i,3],parent_state = node_sim_list[[parent_index]],mu,alpha,site_num,state_num)
  }
  node_character <- t(sapply(node_sim_list,function(x){return(reverse_onehot_encoding(n_state_matrix= x,state_num))}))[1:n_cell,]
  state <- tidyr::unite(data.frame(node_character),col = "state",sep = "")
  character_df <- data.frame(cell = tree$tip.label,state = state$state)
  return(character_df)
}
