# Spatial profiling of colonic immune microenvironment

We performed an unsupervised clustering for immunofluorescence intensity values to determine intensity signatures and to identify each cell according to it. We performed neighboring analysis to see possible cell-cell interactions. 

## Code and data for the application

Required libraries are listed at the beginning of the files. 
Data uploaded in raw_data folder

## Pipeline 
1. Data import and clustering using `WT_non_filtered_outlier_removed.R` code.
2. Neighbouring analysis run using `SPIAT.R` code.
3. Violin plot generation and statistics were performed using `sex violinplots + statistics.R` code.

