## Scripts of Garrido-Oter *et al.*, Modular traits of the Rhizobiales root microbiota and their evolutionary relationship with symbiotic rhizobia, [Cell Host&Microbe, 24(1), 155-167 (2018)](https://www.cell.com/cell-host-microbe/fulltext/S1931-3128(18)30318-4).

originally by Ryohei Thomas Nakano

nakano@mpipz.mpg.de

These scripts are made available to facilitate reproducibility of our research. If you re-use any or part of this code, please reference with comments and cite our paper! Raw data and intermediate results necessary to run these scripts can be downloaded [here](http://www.mpipz.mpg.de/R_scripts).

---------------------------

### Scripts used for processing data, creating the figures and performing the statistical analysis reported in the manuscript.

#### Quality filtering and mapping of sequenced reads:

[1_mapping.sh](https://github.com/rtnakano1984/129E_RNAseq/blob/master/1_mapping.sh)

Bash script for RNAseq read processing and mapping.

[2_makeCountTable.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/2_makeCountTable.R)

R script for preparing read count table from mapped sequence reads.




#### Statistical analysis of gene expression:

[3_GLM_analysis_129E.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/3_GLM_analysis_129E.R)

R script for fitting read counts to a generalized linear model.
[4_DEG_analysis.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/4_DEG_analysis.R)

R script for selecting Differentially Expressed Genes based on Likelihood Ratio Test.
[5_convert_modified_Z.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/5_convert_modified_Z.R)

R script for converting normalized read counts to gene-wise median-centered Z scores.
[6_kmeanClust_AIC.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/6_kmeanClust_AIC.R)

R script for calculating the best number of clusters for k-means.
[7_kmeans.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/7_kmeans.R)

R script for k-mean clustering.
[8_heatmaps.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/8_heatmaps.R)

R script for drawing heatmaps using various values.
[9_enrichment_analysis.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/9_enrichment_analysis.R)

R script for performing enrichment analysis of k-means clusters for Gene Ontology categories.
[10_GLM_analysis_flg22.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/10_GLM_analysis_flg22)

R script for fitting read counts from [Castrillo *et al.*](https://www.nature.com/articles/nature21417).
[11_compare_response-129E_vs_flg22.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/11_compare_response)

R script to compare plant response to *Rhizobium* 129E colonization and a chronic exposure to the flg22 MAMP.
[plotting_parameters.R](https://github.com/rtnakano1984/129E_RNAseq/blob/master/plotting_parameters.R)

R script that contains parameters required for the other scripts.



---------------------------

For any questions regarding these scripts, please contact

Ryohei Thomas Nakano

nakano@mpipz.mpg.de
