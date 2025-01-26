## epigenoproteomics

### Protein Abundance vs RNA and Chromatin Accessibility

A statistical method to identify genes with protein abundance that correlates with RNA abundance and chromatin accessibility (Version: 0.0.0.9000)

Description: 
A statistical method to identify genes with Protein Abundance (LC-MS/MS) that correlate with RNA expression (RNA-seq) and Chromatin Accessibility by ATAC-seq in R. 

It uses a linear model test to assess the probability that the associated chromatin site and RNA expression predict the protein expression 
of the given genes expression and their predicted regulatory elements. 

The first step randomizes the samples into two sets one used for visualization of correlation, and the other used in the linear model. 

The other two functions implement the statistic approach to quantify 
1) the correlation between the pairs of information (RNA-protein and RNA-chromatin) 
2) the protein predictive significance, which is a linear model test.
