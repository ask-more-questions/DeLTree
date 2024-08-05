## DeLTree: A Discrte Model of Lineage Reconstruction with Dynamic Lineage Barcode
Here we introduce a discrete time model to account for the synchronous process of clonal expansion and barcode evolution with a Maximum Likelihood framework. Hence, we model over the cellular geneaology as discrete branching process. For proliferation-active cells like zygote or stem cells, We uses the Galton-Watson process to model over the process of clonal expansion starting from a single progenitor cell after fixed generations,which incorporate the fixed duration of the tracing experiment and provdies intuitive understanding of the edge length as generation time in the reconstructed tree.

### DeLTree Scheme
![Scheme of DeLTree](/SchemeOfDeLTree.jpg)


### Usage 
To reproduce the DeLTree's reconstruction for Subchallenge 1 in DREAM CHALLENGE, we provides the code within the RMarkdown file in DeLTree_on_dream_challenge.html file, and the dataset on dream_challenge_sub1 folder with easy installation of DeLTree package
```
library(devtools)
install_github("ask-more-questions/DeLTree")
library(DeLTree)
```
### Main function and parameters of DeLTree
We provide two function to reconstruct a starting tree from the barcode infomation and a NNI Search function as main functionns for lineage reconstruction. Here we will present the function alongside along with the choice of parameters.The output of tree is versitile for newick and nexus and other acceptable format accepted in `write.tree` in `ape`.

##### Neighbor Joining tree
* `nj_tree(character_info = character_info,site_num = 10,original_state = "1")` <br>
The character_info is a data frame comprising of two variables named cell and barcode. This data frame can be read from format like txt and read.table function in `utils`. The parameter site_number will be the number of sites within one barcode. The paramter "original_state" is default set as "1", which implies the unmutated state on barcode. Here, we wanted to stress that, the temporary function only takes 9 differnt mutant denoted as "0","2",..."9" in one site for the limit of one hot encoding.

##### Cherry Reconstruction
* `likelihood_based_recon(character_info,mu,alpha,non_bifur_pro,alternative_threshold,nGen,site_num,state_num,precluster_replicate)` <br>
This above function describes a reconstruction method with Cherry as basic unit. In which we first reconstruct a the pairwise likelihood ratio of having one or none bifurcation in two generation time. The implied barcode of Most Recent Common Ancentor vertex follows the irreversible parsimony rule. Here, `mu` and `alpha` will be the site-specific mutation rate and editing priors on one generation time in barcode. In our model, the growth characteristic of lineages lies in the parameter,`non_bifur_pro`, which will be the proportion of cell not bifurcated after one generation time. Alternative_threshold will be the set default as zero and if 1,two sequences are 10 times more likely to have a bifurcation event in 2 generation. Besides, the prameter `precluster_replicate` is a logical parameter indicating whether or not to cluster the cells with identical barcodes into balanced cluster before reconstruct phylogeny.


##### DeLTree Search with Discrete Edge Length
* `nni_iter_withedgelength_pseudonode(current_tree,mu,alpha,nGen,non_bifur_pro,state_num,edgelength_assignment)` <br>
In this function, we take  a tree function of phylo structure with its tiplabel formatted as 'cell_barcode' as input for current tree. Here the `nGen` implies a prefixed tree height based on experimental duration on cell divisions. The paramter `edgelength_assignment` provides two options as "bottom up iteration" or "direct assignment", where the former applies the longest pending edge length as initiation and performs bottom-up approach to local optimization, while the latter utilizes the depth of node and discrete edge length rule as constraints for edge length of cherry structure, which is used as default.

