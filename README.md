# violinEnsemble
Create violin plots of gene expression for multiple genes

Given a Seurat object and a list of genes, groupedViolins will plot expression levels for each gene with the cells 
grouped by a provided variable (such as cell type or cluster identity).  If the gene list provided is a table containing the
columns "cluster" and "gene" (such as from FindAllClusters(), the genes will be grouped by the cluster with which 
they are associated.

As an example, here are the top 3 markers for each cluster in the dataset from [Villani et al. Science 2017 paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5775029/)
![Villani et.al. plot](https://github.com/milescsmith/violinEnsemble/blob/master/example.jpeg)
