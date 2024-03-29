% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_metric.R
\name{triplet_distance}
\alias{triplet_distance}
\title{Tree Metric}
\usage{
triplet_distance(tree1, tree2, normalized = TRUE)
}
\arguments{
\item{tree1}{An oject of class "phylo"}

\item{tree2}{An oject of class "phylo"}

\item{normalized}{default set to TRUE; to normalize the triplet distance by the total number of possible triplets.}
}
\value{
Triplet distance between the two trees
}
\description{
Function to calculate triplet distance between two phylogenetic trees. The depths of MRCA and other relevant nodes are used to identify the cherry in each tree.
It takes two phylogenetic trees andcalculates the triplet distance between two phylogenetic trees.
The triplet distance measures the number of triplets (unordered sets of three taxa)
that differ in topology between the two trees.
}
