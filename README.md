## DeLTree: A Discrte Model of Lineage Reconstruction with Dynamic Lineage Barcode
Here we introduce a discrete time model to account for the synchronous process of clonal expansion and barcode evolution with a Maximum Likelihood framework. Hence, we model over the cellular geneaology as discrete branching process. For proliferation-active cells like zygote or stem cells, We uses the Galton-Watson process to model over the process of clonal expansion starting from a single progenitor cell after fixed generations,which incorporate the fixed duration of the tracing experiment and provdies intuitive understanding of the edge length as generation time in the reconstructed tree.

### DeLTree Scheme
[![Scheme of DeLTree]](www.github.com/ask-more-questions/DeLTree/main/SchemeOfDeLTree.jpg)


### Usage 
To reproduce the DeLTree's reconstruction for Subchallenge 1 in DREAM CHALLENGE, we provides the code within the RMarkdown file in DeLTree_on_dream_challenge.html file, and the dataset on dream_challenge_sub1 folder with easy installation of DeLTree package
```
library(devtools)
install_github("ask-more-questions/DeLTree")
library(DeLTree)
```

