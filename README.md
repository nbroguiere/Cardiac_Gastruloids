# Cardiac_Gastruloids
Single-cell RNA-seq analysis of cardiac gastruloids for [Rossi et al. Cell Stem Cell 2020](https://www.cell.com/cell-stem-cell/pdf/S1934-5909(20)30507-5.pdf).

# Data exploration (recommended)
Download the data from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158999), install R-studio and [Seurat](https://satijalab.org/seurat/) v3, and use the file Explore_Data.R for a quick dive into cell annotations and gene expression, and to be setup for further exploration in R. 

# Complete analysis reproduction (advanced users)
To reproduce the complete analysis, obtain the raw data from Rossi et al. 2020 (cardiac gastruloids, fastq files on SRA) as well as Pijuan-Sala et al. 2019 (used as a reference in vivo atlas). Install all dependencies listed in SessionInfo.txt, download the input tables, cell-type classifier, and complete code, and run cellranger for alignment, [velocyto](http://velocyto.org/) for RNA splicing estimates, followed by CardiacGastruloidsAnalysis.R for the main analysis, and by the jupyter lab notebooks for RNA-velocity calculations with [scVelo](https://github.com/theislab/scvelo). 
