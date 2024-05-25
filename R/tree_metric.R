#' Tree Metric
#'
#' @param tree1  An oject of class "phylo"
#' @param tree2  An oject of class "phylo"
#' @param normalized  default set to TRUE; to normalize the triplet distance by the total number of possible triplets.
#'
#' @return Triplet distance between the two trees
#'
#' @description
#' Function to calculate triplet distance between two phylogenetic trees. The depths of MRCA and other relevant nodes are used to identify the cherry in each tree.
#'     It takes two phylogenetic trees andcalculates the triplet distance between two phylogenetic trees.
#'     The triplet distance measures the number of triplets (unordered sets of three taxa)
#'     that differ in topology between the two trees.
#' @export
triplet_distance <- function(tree1, tree2,normalized = TRUE) {
  # Get the number of taxa in the trees
  n_taxa <- length(tree1$tip.label)
  triplet_dist <- 0
  # Calculate node depths for each tree
  nodedepth1 <- node.depth(tree1,method = 1)
  nodedepth2 <- node.depth(tree2,method = 1)
  # Iterate through all possible triplets in the trees
  for (i in 1:(n_taxa - 2)) {
    for (j in (i+1):(n_taxa - 1)) {
      for (k in (j+1):n_taxa) {
        # Extract and sort the triplet labels
        triplet <- sort(c(tree1$tip.label[i], tree1$tip.label[j], tree1$tip.label[k]),decreasing = TRUE)
        # Get depths of MRCA (Most Recent Common Ancestor) and other relevant nodes for tree1
        triplet1_mrca_depth<- nodedepth1[getMRCA(tree1,tip= triplet)]
        mrca1_depth <- nodedepth1[getMRCA(tree1,tip = triplet[c(1,2)])]
        mrca2_depth <- nodedepth1[getMRCA(tree1,tip = triplet[c(2,3)])]
        mrca3_depth <- nodedepth1[getMRCA(tree1,tip = triplet[c(1,3)])]
        # Identify the cherry (pair of descendants) in tree1
        if (mrca1_depth != triplet1_mrca_depth){
          triplet1 <- c(triplet[3],triplet[c(1,2)])
        }else if  (mrca2_depth != triplet1_mrca_depth){
          triplet1 <- c(triplet[1],triplet[c(2,3)])
        }else if (mrca3_depth != triplet1_mrca_depth){
          triplet1 <- c(triplet[2],triplet[c(1,3)])
        }else print("Error: cannot find the cherry")
        # same process for tree two
        triplet2_mrca_depth <- nodedepth2[getMRCA(tree2,tip= triplet)]
        mrca1_depth <- nodedepth2[getMRCA(tree2,tip = triplet[c(1,2)])]
        mrca2_depth <- nodedepth2[getMRCA(tree2,tip = triplet[c(2,3)])]
        mrca3_depth <- nodedepth2[getMRCA(tree2,tip = triplet[c(1,3)])]
        if (mrca1_depth != triplet2_mrca_depth){
          triplet2 <- c(triplet[3],triplet[c(1,2)])
        }else if  (mrca2_depth != triplet2_mrca_depth){
          triplet2 <- c(triplet[1],triplet[c(2,3)])
        }else if (mrca3_depth != triplet2_mrca_depth){
          triplet2 <- c(triplet[2],triplet[c(1,3)])
        }else print("Error: cannot find the cherry")
        # Check if the triplets have the same topology in both trees
        if(!identical(triplet1,triplet2))
          triplet_dist <- triplet_dist + 1
      }
    }
  }
  if (normalized)
    triplet_dist <- triplet_dist/choose(n_taxa,3)
  return(triplet_dist)
}


#' Weighted Tree Metric
#'
#' @param tree1  An oject of class "phylo"
#' @param tree2  An oject of class "phylo"
#' @param normalized default set to TRUE; to normalize the triplet distance by the total number of possible triplets.
#'
#' @description
#'     This function extends the triplet_distance function by incorporating information about the branch lengths of the phylogenetic trees.
#'     For matching triplets(identical topology),the weight is calculated as the absolute difference in edge lengths between the common ancestor of the triplet in tree1 and tree2.
#'     For Non-matching Triplets,the weight is calculated as the average edge length of the common ancestor in both trees.
#' @return triplet distance between the two trees
#'
#' @export

triplet_distance_weighted <- function(tree1, tree2,normalized = TRUE) {
  # Get the number of taxa in the trees
  n_taxa <- length(tree1$tip.label)
  # Initialize weighted triplet distance
  triplet_dist <- 0
  # Calculate node depths and edge lengths for each tree
  nodedepth1 <- node.depth(tree1,method = 1)
  nodedepth2 <- node.depth(tree2,method = 1)
  nodedepth_edgelength1 <- node.depth.edgelength(tree1)
  nodedepth_edgelength2 <- node.depth.edgelength(tree2)
  # Iterate through all possible triplets in the trees
  for (i in 1:(n_taxa - 2)) {
    for (j in (i+1):(n_taxa - 1)) {
      for (k in (j+1):n_taxa) {
        # Extract and sort the triplet labels
        triplet <- sort(c(tree1$tip.label[i], tree1$tip.label[j], tree1$tip.label[k]),decreasing = TRUE)

        # Get depths of MRCA (Most Recent Common Ancestor) and other relevant nodes for tree1
        triplet1_mrca_depth<- nodedepth1[getMRCA(tree1,tip= triplet)]
        mrca1_depth <- nodedepth1[getMRCA(tree1,tip = triplet[c(1,2)])]
        mrca2_depth <- nodedepth1[getMRCA(tree1,tip = triplet[c(2,3)])]
        mrca3_depth <- nodedepth1[getMRCA(tree1,tip = triplet[c(1,3)])]
        # Identify the cherry (pair of descendants) in tree1
        if (mrca1_depth != triplet1_mrca_depth){
          triplet1 <- c(triplet[3],triplet[c(1,2)])
        }else if  (mrca2_depth != triplet1_mrca_depth){
          triplet1 <- c(triplet[1],triplet[c(2,3)])
        }else if (mrca3_depth != triplet1_mrca_depth){
          triplet1 <- c(triplet[2],triplet[c(1,3)])
        }else print("Error: cannot find the cherry")

        # Get depths of MRCA and other relevant nodes for tree2
        triplet2_mrca_depth <- nodedepth2[getMRCA(tree2,tip= triplet)]
        mrca1_depth <- nodedepth2[getMRCA(tree2,tip = triplet[c(1,2)])]
        mrca2_depth <- nodedepth2[getMRCA(tree2,tip = triplet[c(2,3)])]
        mrca3_depth <- nodedepth2[getMRCA(tree2,tip = triplet[c(1,3)])]

        # Check if the triplets have the same topology in both trees
        if (mrca1_depth != triplet2_mrca_depth){
          triplet2 <- c(triplet[3],triplet[c(1,2)])
        }else if  (mrca2_depth != triplet2_mrca_depth){
          triplet2 <- c(triplet[1],triplet[c(2,3)])
        }else if (mrca3_depth != triplet2_mrca_depth){
          triplet2 <- c(triplet[2],triplet[c(1,3)])
        }else print("Error: cannot find the cherry")
        if (identical(triplet1,triplet2)){
          weight_tree1as_reference <- abs(nodedepth_edgelength1[getMRCA(tree1,tip= triplet)] - nodedepth_edgelength2[getMRCA(tree2,tip= triplet)])
          triplet_dist <- triplet_dist + 1*weight_tree1as_reference
        } else{
          weight_tree1as_reference <- (abs(nodedepth_edgelength1[getMRCA(tree1,tip= triplet)] + nodedepth_edgelength2[getMRCA(tree2,tip= triplet)]))/2
          triplet_dist <- triplet_dist + 1*weight_tree1as_reference
        }
      }
    }
  }
  # Normalize weighted triplet distance
  if (normalized)
    triplet_dist <- triplet_dist/choose(n_taxa,3)
  return(triplet_dist)
}


#' Compute Colless imbalance
#'
#' @param tree a phylogenetic tree
#'
#' @return colless imbalance value

colless_index  <- function(tree){
  L <- length(tree$tip.label)
  colless_value <- 0
  for(i in 1:tree$Nnode){
    child_taxa <- Descendants(x=tree,node = L+i,type = "children")
    if(length(child_taxa) == 1)
      next
    differ <- abs(node.depth(phy = tree,method = 1)[child_taxa[1]]-node.depth(phy = tree,method = 1)[child_taxa[2]])
    colless_value <- colless_value+differ
  }
  colless_value <- colless_value*(2/((L-1)*(L-2)))
  return(colless_value)
}

#' Weighted colless index
#'
#' @param tree a phylogenetic tree
#'
#' @return  colless imbalance value with the weight descirbed in the reference paper

colless_index_weighted <- function(tree){
  L <- length(tree$tip.label)
  colless_value <- 0
  for(i in 1:tree$Nnode){
    child_taxa <- Descendants(x=tree,node = L+i,type = "children")
    if(length(child_taxa) == 1)
      next
    differ <- abs(node.depth(phy = tree,method = 1)[child_taxa[1]]-node.depth(phy = tree,method = 1)[child_taxa[2]])
    g <-  max(node.depth(tree,method = 2))-node.depth(tree,method = 2)[L+i]
    colless_value <- colless_value+differ*(1/2^g)
  }
  N <- 0
  for(j in 0:(L-3)){
    N <- N+(L-2-j)/2^j
  }
  colless_w <- colless_value/N
  return(colless_w)
}
