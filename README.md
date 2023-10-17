# AD-TransComp-R
# The following scripts were used for implementation of TransComp-R with AD single-cell RNA-sequencing Data
# Step 1: Ortholog Finder.Rmd - finds orthologs between mouse and human data using orthogene package
# Step 2: AD Mouse Processing Orthologs - processes mouse AD data
# Step 3: AD_human_rscript - using hpc, run initial processing of human AD data due to large size
# Step 4: AD Mouse Processing Orthologs - finish human AD data processing
# Step 5: AD Ortholog PCA - processes human and mouse microglia so the cell types are accurate, runs PCA implementation of TransComp-R
# Step 6: AD Ortholog sPCA - runs sPCA implementation of TransComp-R
# Step 7: AD Differential Expression Microglia - differential expression for mouse and human microglial cells
# Step 8: AD GSEA RAnked - fGSEA analysis on significant PCs from PCA and sPCA analysis
# Step 9: AD Visualization - script for figure identification
