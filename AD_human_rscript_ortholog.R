# Library Import
library(BiocManager, lib.loc = "/home/axb1657/R-4.1.1")
library(Seurat)
library(tidyverse, lib.loc = "/home/axb1657/R-4.1.1")
library(patchwork) 
library(dittoSeq, lib.loc = "/home/axb1657/R-4.1.1") 
library(stringr)
library(gprofiler2, lib.loc = "/home/axb1657/R-4.1.1")
library(harmony)
library(patchwork)
library(cowplot)
library(openxlsx)
library(glmGamPoi, lib.loc = "/home/axb1657/R-4.1.1")
#library(MAST, lib.loc = "/home/axb1657/R-4.1.1")


home_dir <- "/home/axb1657/Brubaker_Summer_Research_2023"


# Data Import

samples <- c("A1", "A2", "A3", "A4")
names <- c("AD_1", "C_1", "AD_2", "C_2")
data_list <- list()

# Load in human data
for(i in samples){ 
  data <- Read10X(paste0(home_dir, "/AD_human_raw/", i)) 
  data_list[[i]] <- data
} 

human_orthologs <- readRDS(file = paste0(home_dir, "/AD_human_objects_ortholog/human_orthologs.RData"))

# Create Seurat object of human orthologs
for(i in samples){ 
  data_list[[i]] <- CreateSeuratObject(counts = data_list[[i]], project = i)
} 

# Filter Seurat objects based on human orthologs
for(i in samples){
  data_list[[i]] <- data_list[[i]][human_orthologs]
}

names(data_list) <- names


# Quality Control

for(i in names){ 
  data_list[[i]]@meta.data$Sample <- i
}

# Mitochondrial quality control
for(i in names){ 
  data_list[[i]][["percent.mt"]] <- PercentageFeatureSet(data_list[[i]], pattern = "^MT-") 
} 

# RNA quality control
for(i in names){ 
  data_list[[i]][["rna.size"]] <-   log10(Matrix::colSums(data_list[[i]]@assays$RNA@counts)) 
} 


#QC Visualization 

combined_test <- merge(x = data_list[[1]], y = data_list[2:length(data_list)], merge.data=TRUE)
p1 <- VlnPlot(combined_test, "rna.size", group.by = "Sample")
#print(p1)

p2 <- VlnPlot(combined_test, "percent.mt", group.by = "Sample")
#print(p2)

cairo_ps(filename= paste0(home_dir, "/Human_Viz_ortholog/S1C.eps"), width = 7, height = 7) # name the file as an EPS
p1 #variable of the figure that is defined/saved, can be used any name
dev.off()

cairo_ps(filename= paste0(home_dir, "/Human_Viz_ortholog/S1D.eps"), width = 7, height = 7) # name the file as an EPS
p2 #variable of the figure that is defined/saved, can be used any name
dev.off()

rm(combined_test)

for(i in names){ 
  qc_plot <- VlnPlot(data_list[[i]], features = c("rna.size", "percent.mt"), ncol = 3) 
  print(qc_plot)
  ggsave(paste0(home_dir, "/Human_Viz_ortholog/", i, "_qc.pdf"))
} 


#10% for mitochondrial filtering threshold
#+/- 3 median absolute deviations from the median rna library size



rna_size_min <- list() 
rna_size_max <- list() 


# Calculates min and max for rna filtering
for(i in names){ 
  rna_size_min[[i]] <- median(data_list[[i]]@meta.data$rna.size) - (3*mad(data_list[[i]]@meta.data$rna.size)) 
  
  rna_size_max[[i]] <- median(data_list[[i]]@meta.data$rna.size) + (3*mad(data_list[[i]]@meta.data$rna.size)) 
} 


# Save QC metric DF
qc_metrics <- c() 
for(i in names){ 
  preqc_cells <- ncol(data_list[[i]]) 
  qc_metrics <- rbind(qc_metrics, preqc_cells) 
}

colnames(qc_metrics) <- "preqc_cells" 
row.names(qc_metrics) <- NULL 
preqc_metrics <- as.data.frame(qc_metrics)


# Filter based on QC metrics
for(i in names){ 
  temp_max <- rna_size_max[[i]] 
  temp_min <-  rna_size_min[[i]] 
  data_list[[i]] <- subset(data_list[[i]], subset = percent.mt < 10 & rna.size < temp_max & rna.size > temp_min) 
  
} 


# Finish QC Table
qc_metrics <- c() 
for(i in names){ 
  preqc_cells <- ncol(data_list[[i]]) 
  qc_metrics <- rbind(qc_metrics, preqc_cells) 
} 



colnames(qc_metrics) <- "postqc_cells" 
row.names(qc_metrics) <- NULL 
postqc_metrics <- as.data.frame(qc_metrics)
cellcounts <- cbind(preqc_metrics, postqc_metrics)
write.csv(cellcounts, file = paste0(home_dir, "/human_cellcounts_ortholog.csv"))


# Cell Cycle Scoring Addition

s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes 

for(i in names){ 
  data_list[[i]] <- NormalizeData(data_list[[i]]) 
  data_list[[i]] <- CellCycleScoring(data_list[[i]], s.features = s.genes, g2m.features = g2m.genes) 
} 

data_list$AD_1@meta.data$Genotype <- "AD"
data_list$AD_2@meta.data$Genotype <- "AD"
data_list$C_1@meta.data$Genotype <- "C"
data_list$C_2@meta.data$Genotype <- "C"


# SCTransform

#Applies SCTransform individually to each sample 

data_list <- lapply(X = data_list, FUN = SCTransform, method = "glmGamPoi", return.only.var.genes = FALSE, vars.to.regress = c("S.Score", "G2M.Score")) 




#saveRDS(data_list, file = paste0(home_dir, "/AD_human_objects/combined.RData"))
#data_list <- readRDS(file = paste0(home_dir, "/AD_human_objects/combined.RData"))


# Selects variable features across samples
variable.features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000) 



# Merge, batch correction
data_list <- merge(x = data_list[[1]], y = data_list[2:length(data_list)], merge.data=TRUE) 

VariableFeatures(data_list) <- variable.features


# Run PCA
data_list <- RunPCA(data_list, verbose = FALSE) 

# Variances before integration
stdev <- data_list@reductions$pca@stdev
vaiances <- (stdev)^2 / sum(stdev^2)
write.csv(variances, file = paste0(home_dir, "/pca_variances.csv"))



#Checking differences in principal components 

#options(repr.plot.height = 12, repr.plot.width = 20) 
p1 <- DimPlot(object = data_list, reduction = "pca", pt.size = .1, group.by = "Sample") 
#p2 <- VlnPlot(object = data_list, features = "PC_1", group.by = "orig.ident", pt.size = 0) 
#plot_grid(p1,p2) 
cairo_ps(filename= paste0(home_dir, "/Human_Viz_ortholog/1F_pre.eps"), width = 9, height = 7) # name the file as an EPS
p1 #variable of the figure that is defined/saved, can be used any name
dev.off()
#ggsave(paste0(home_dir, "/Human_Viz_ortholog/pca_prebatch.pdf"))

#Batch correction necessary in this instance

#Running Harmony  
#options(repr.plot.height = 8, repr.plot.width = 16) 

data_list <- RunHarmony(data_list, assay.use = "SCT", group.by.vars = "Sample", plot_convergence = TRUE) 
#ggsave(paste0(home_dir, "/Human_Viz_ortholog/convergence_plot.pdf"))


#Visualizing batch corrected data 

p1 <- DimPlot(object = data_list, reduction = "harmony", pt.size = .1, group.by = "Sample") 

#p2 <- VlnPlot(object = data_list, features = "harmony_1", group.by = "orig.ident", pt.size = 0) 

#plot_grid(p1,p2) 
cairo_ps(filename= paste0(home_dir, "/Human_Viz_ortholog/1F_post.eps"), width = 9, height = 7) # name the file as an EPS
p1 #variable of the figure that is defined/saved, can be used any name
dev.off()
#ggsave(paste0(home_dir, "/Human_Viz_ortholog/pca_postbatch.pdf"))

stdev <- data_list@reductions$harmony@stdev
variances <- (stdev)^2 / sum(stdev^2)
write.csv(variances, file = paste0(home_dir, "/harmony_variances.csv"))

#Looking better.

ElbowPlot(data_list, reduction = "harmony", n = 50) 
ggsave(paste0(home_dir, "/Human_Viz_ortholog/elbowplot.pdf"))


#Lets try 30 PCs

# Reducing Dimensionality
#Performing Dimensional Reduction 

data_list <- RunUMAP(data_list, reduction = "harmony", dims = 1:30, verbose = FALSE) 

data_list <- FindNeighbors(data_list, reduction = "harmony", dims = 1:30, verbose = FALSE) 

data_list <- FindClusters(data_list, resolution = 0.6, verbose = FALSE) 


# Visualization
dittoDimPlot(data_list, "seurat_clusters", do.label = TRUE, legend.show = FALSE, labels.size = 3) 
ggsave(paste0(home_dir, "/Human_Viz_ortholog/seurat_clusters.pdf"))

dittoDimPlot(data_list, "Genotype") 
ggsave(paste0(home_dir, "/Human_Viz_ortholog/genotype.pdf"))

dittoDimPlot(data_list, "orig.ident") 
ggsave(paste0(home_dir, "/Human_Viz_ortholog/orig_ident.pdf"))

dittoDimPlot(data_list, "percent.mt") 
ggsave(paste0(home_dir, "/Human_Viz_ortholog/percent_mt.pdf"))

dittoDimPlot(data_list, "Phase") 
ggsave(paste0(home_dir, "/Human_Viz_ortholog/phase.pdf"))



# Save object
saveRDS(data_list, file = paste0(home_dir, "/AD_human_objects_ortholog/human_adc_final.RData"))
