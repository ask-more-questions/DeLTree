## DeLTree: A Discrte Model of Lineage Reconstruction with Dynamic Lineage Barcode
Here we introduce a discrete time model to account for the synchronous process of clonal expansion and barcode evolution with a Maximum Likelihood framework. Hence, we model over the cellular geneaology as discrete branching process. For proliferation-active cells like zygote or stem cells, We uses the Galton-Watson process to model over the process of clonal expansion starting from a single progenitor cell after fixed generations,which incorporate the fixed duration of the tracing experiment and provdies intuitive understanding of the edge length as generation time in the reconstructed tree.

### DeLTree Scheme
![Scheme of DeLTree](/SchemeOfDeLTree.jpg)


### Usage 
To reproduce the DeLTree's reconstruction for Subchallenge 1 in DREAM CHALLENGE, we provides the code within the RMarkdown file in DeLTree_scripts.R file, and the dataset on dream_challenge_sub1 folder with easy installation of DeLTree package
```
library(devtools)
install_github("ask-more-questions/DeLTree")
library(DeLTree)
library(caret)
```
### Main function and parameters of DeLTree
We provide two function to reconstruct a starting tree from the barcode infomation and a NNI Search function as main functionns for lineage reconstruction. Here we will present the function alongside along with the choice of parameters.The output of tree is versitile for newick and nexus and other acceptable format accepted in `write.tree` in `ape`.

##### Neighbor Joining tree
* `nj_tree(character_info = character_info,site_num = 10,original_state = "1")` <br>
The character_info is a data frame comprising of two variables named cell and barcode.This data frame can be read from format like txt and read.table function in `utils`. <br>
The parameter `site_number` is the number of sites within one barcode. The paramter "original_state" is default set as "1", which implies the unmutated state on barcode. Here, we wanted to stress that, the temporary function only takes 9 differnt mutant denoted as "0","2",..."9" in one site for the limit of one-hot encoding step we use in `caret` package.

##### Cherry Reconstruction
* `cherry_based_recon(character_info,mu,alpha,non_bifur_pro,alternative_threshold,nGen,site_num,state_num,precluster_replicate)` <br>
This above function describes a reconstruction method with Cherry as basic unit. <br>
Here, `mu` and `alpha` are the site-specific mutation rate and editing priors on one generation time in barcode respectively. 
`non_bifur_pro` is the proportion of cell not bifurcated after one generation time, reepresenting growth characteristic of lineages. `Alternative_threshold` is the set default as zero and if 1,two sequences are 10 times more likely to have a bifurcation event in 2 generation.


##### DeLTree Search with Discrete Edge Length
* `nni_deltree(current_tree,mu,alpha,nGen,non_bifur_pro,state_num,edgelength_assignment)` <br>
In this function, we take  a tree function of phylo structure with its tip label formatted as 'cell_barcode' as input for current tree. and perform one nni move and output the best tree as a phylo structure and the corresponding likelihood. <br>
Here the `nGen` implies a prefixed tree height based on experimental duration scaled on cell divisions. Paramter `edgelength_assignment` provides two options as "bottom up iteration" or "direct assignment", where the former applies the longest pending edge length as initiation and performs bottom-up approach to local optimization, while the latter utilizes the depth of node and discrete edge length rule as constraints for edge length of cherry structure, which is used as default.

### Run
#### Discrete bifurcation on tree simulation
Here we present a simulation example with fixed time duration and bifurcation probability per generation time.<br>
Supposing a no death scheme for $g$ generation time and expected clonal size of $N$ cells. We compute an estimated bifurctaion probability per generation time as $p_{bifur} = \exp(\frac{log(N)}{g})-1$. Here we estimate the $p_{bifur}$ as 0.76 when $g$ = 6, $N$ = 30. <br> 
```
sim_tree = topology_simulation(bifur_pro = 0.76 ,nGen = 6)
mu = rep(0.15,10)
alpha = rep(list(0.5,0.5),10)
lineage = lineage_sim(tree = sim_tree,state_num = 3, site_num = 10,mu,alpha)
sim_tree$tip.label <- sapply(sim_tree$tip.label,function(x){
    state_index <-lineage$cell %in% x
    cell_state <- paste(x,lineage$state[state_index],sep = "_")
    return(cell_state)})
```
The below Figure shows an example of a simulated tree with the setting above.

![image](https://github.com/user-attachments/assets/fb59f8fe-dd48-4436-908d-173185e1b540)

#### Run Neighbor Joining and DeLTree NNI
For barcode information saved in sub_test_1.txt in the dream_challenge_sub1 folder.
```
lineage_file = "dream_challenge_sub1/sub1_test_1.txt"
character_info = read.table(file = lineage_file,header= TRUE,colClass = "character")
NJ_tree = nj_tree(character_info = character_info,site_num = 10,original_state = "1")
max_iter = 10
mu = rep(0.15,10)
alpha = rep(list(0.5,0.5),10)
non_bifur_pro = 0.2
current_tree_del <- direct_assignment(phylo = NJ_tree,nGen = max(node.depth(NJ_tree,method = 2)),state_num = 3,mu = mu,alpha = alpha,non_bifur_pro = non_bifur_pro)
current_likelihood <- deltree_likelihood(discrete_EdgeLength_tree = current_tree_del,state_num = 3, mu = mu,alpha = alpha,non_bifur_pro = non_bifur_pro)[3]
for (j in (1:max_iter)){
  nni_info <- nni_deltree(current_tree = current_tree_del,mu,alpha,nGen=max(node.depth(current_tree_del,method = 2)),non_bifur_pro,state_num =3,edgelength_assignment = "direct assignment")
  if (max(nni_info$likelihood)>current_likelihood){
    current_tree_del <- nni_info$best_tree
    current_likelihood <- max(nni_info$likelihood)
  } else break
}
NJ_DelTree = current_tree_del
```
Here parameter `nGen` is set as as the topological height of the current tree as a loosen tree height range. It is recommended as a default setting to coincide with the non-zero discrete edge length assumption. 
