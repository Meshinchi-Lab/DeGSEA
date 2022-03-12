# R package for Differential Expression and GSEA using RNA-seq
### for DeGSEA

### Installation 

Install using `devtools` 
```
devtools::install_github("Meshinchi-Lab/DeGSEA")
```

The package was developed for use in the Meshinchi lab to help streamline the association of clinical covariates and RNA-seq expression data. The clincal data elements are can be found at the [TARGET data matrix](https://ocg.cancer.gov/programs/target/data-matrix). The RNA-sequencing data can be found at the [Genomic Data Commons](https://portal.gdc.cancer.gov/) under project ID "TARGET-AML". 

Its aim is only to simplify the workflows of differential expression analysis, GSEA, and unsupervised clustering methods, as well the data vizualization using open source, publicly avaialble R packages such as `limma` and `edgeR`. 

Note that some functions are fairly outdated, but are retained for backwards compatibility of analyses run from ~2017 - 2018. 

Author: Jenny Leopoldina Smith<br>
ORCID: [0000-0003-0402-2779](https://orcid.org/0000-0003-0402-2779)
<br>
