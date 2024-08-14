library(devtools)
install_github("ask-more-questions/DeLTree")
library(DeLTree)
library(phangorn)
library(tidyr)
library(caret)
filepath_train_barcode <- paste("dream_challenge_sub1/groundtruth_train/","sub1_train_",1:76,".txt",sep = "")
filepath_train_tree <- paste("dream_challenge_sub1/groundtruth_train/","sub1_train_",1:76,".nwk",sep = "")
filepath_test_barcode <- paste("dream_challenge_sub1/","sub1_test_",1:30,".txt",sep = "")
GroundTruth_tree <- read.tree("dream_challenge_sub1/Goldstandard_SC1.txt")



# Estimation of two parameter : non_bifurcation_pro and mutation rate with training dataset
site_num <- 10
state_num <- 3
branch_1 <- c()
branch_2p <- c()
for(i in 1:length(filepath_train_tree)){
  train_tree <- read.tree(filepath_train_tree[i])
  train_tree$edge.length[train_tree$edge[,2] <= length(train_tree$tip.label)] <- 0
  branch_1[i] <-sum(round(train_tree$edge.length/48) >= 1)
  branch_2p[i] <- sum(round(train_tree$edge.length/48)[round(train_tree$edge.length/48) >=2]-1)
}
non_bifur_pro <- 1-sum(branch_1)/sum(branch_1+branch_2p)

count_info <- matrix(data = 0,nrow = state_num,ncol = site_num)
for (i in 1:length(filepath_train_tree)){
  tree <- read.tree(file = filepath_train_tree[i])
  tree$edge.length <- round(tree$edge.length/48)
  count_info <- count_info + get_site_mutation(tree,site_num,state_num)
}

mu <- (apply(count_info,2,sum)-count_info[2,])/apply(count_info,2,sum)
alpha <- list()
for(i in 1:site_num){
  alpha[[i]] <- sapply(1:state_num,function(x){return(count_info[x,]/(apply(count_info,2,sum)-count_info[2,]))})[,-2][i,]
  if(round(sum(alpha[[i]]),digits = 2) !=1)
    stop("sum of alpha[[i]] doesnot equal to 1 !")
}

## Starting tree reconstruction: Neighbor Joining
start_topo <- list()
start_tree_del <- list()
for (i in 1:length(filepath_test_barcode)){
  character_info <- read.table(file = filepath_test_barcode[i],header = TRUE,colClasses = "character")
  start_topo <- nj_tree(character_info,site_num,original_state = "1")
  #cheery reconstruction
  #start_topo <- cherry_based_recon(character_info,mu,alpha,non_bifur_pro,alternative_threshold = 0,site_num,state_num,precluster_replicate =TRUE)
  start_topo_del <- direct_assignment(phylo = start_topo,nGen =max(node.depth(start_topo,method = 2)),state_num =3,mu,alpha,non_bifur_pro)
  start_tree_del[[i]] <- start_topo_del
}


recon_rf_d <- c()
recon_tri_d <- c()
for (i in 1:30){
  recon_rf_d[i] <- RF.dist(tree1 =GroundTruth_tree[[i]],tree2 = start_tree_del[[i]],normalize = TRUE,rooted = TRUE )
  recon_tri_d[i] <-triplet_distance(tree1 =GroundTruth_tree[[i]],tree2 = start_tree_del[[i]],normalized = TRUE)
}
print(colMeans(cbind(recon_rf_d,recon_tri_d)))

## NNI Tree Search : DeLTree NNI

current_tree <- edgelength_assignment <-"direct assignment"
nni_recorder <- list()
for (i in 1:30){
  current_tree_del <- start_tree_del[[i]]
  current_likelihood <- deltree_likelihood(discrete_EdgeLength_tree = current_tree_del,state_num = 3,mu = mu,alpha = alpha,non_bifur_pro = non_bifur_pro)[3]
  move <- 0
  for (j in 1:10){
    nni_info <- nni_deltree(current_tree_del,mu,alpha,nGen=max(node.depth(remove_pseudonode(current_tree_del),method = 2)),non_bifur_pro,state_num =3,edgelength_assignment = edgelength_assignment)
    if (max(nni_info$likelihood)>current_likelihood){
      current_tree_del <- nni_info$best_tree
      current_likelihood <- max(nni_info$likelihood)
      move <- move + 1
    } else break
  }
  nni_recorder$best_score[i] <- current_likelihood
  nni_recorder$best_tree[[i]] <- current_tree_del
  nni_recorder$move[i] <- move
}


recon_nni_rf_d <- c()
recon_nni_tri_d <- c()
for (i in 1:30){
  recon_nni_rf_d[i] <- RF.dist(tree1 =GroundTruth_tree[[i]],tree2 = nni_recorder$best_tree[[i]],normalize = TRUE,rooted = TRUE )
  recon_nni_tri_d[i] <-triplet_distance(tree1 =GroundTruth_tree[[i]],tree2 = nni_recorder$best_tree[[i]],normalized = TRUE)
}
print(colMeans(cbind(recon_nni_rf_d,recon_nni_tri_d)))

