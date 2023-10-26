# AD-TransComp-R
## The following scripts were used for implementation of TransComp-R with AD single-cell RNA-sequencing Data <br>
Step 1: Ortholog Finder.Rmd - finds orthologs between mouse and human data using orthogene package <br>
Step 2: AD Mouse Processing Orthologs - processes mouse AD data <br>
Step 3: AD_human_rscript_ortholog - using hpc, run initial processing of human AD data due to large size <br>
Step 4: AD Human Processing Orthologs - finish human AD data processing <br>
Step 5: AD Ortholog PCA - processes human and mouse microglia so the cell types are accurate, runs PCA implementation of TransComp-R <br>
Step 6: AD Ortholog sPCA - runs sPCA implementation of TransComp-R <br>
Step 7: Therapeutic Testing of principal components <br>

Other scripts <br>
AD GSEA Ranked: GSEA analysis for signifcant PCs in PCA and sPCA <br>
AD Differential Expression Microglia: DE for Microglia

