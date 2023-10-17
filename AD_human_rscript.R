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
names <- c("AD1", "C1", "AD2", "C2")
data_list <- list()

for(i in samples){ 
  data <- Read10X(paste0(home_dir, "/AD_human_raw/", i)) 
  data_list[[i]] <- data 
} 


for(i in samples){ 
  data_list[[i]] <- CreateSeuratObject(counts = data_list[[i]], project = i) 
} 
names(data_list) <- names



# Quality Control

for(i in names){ 
  data_list[[i]][["percent.mt"]] <- PercentageFeatureSet(data_list[[i]], pattern = "^MT-") 
} 

for(i in names){ 
  data_list[[i]][["rna.size"]] <-   log10(Matrix::colSums(data_list[[i]]@assays$RNA@counts)) 
} 


#QC Visualization 

for(i in names){ 
  qc_plot <- VlnPlot(data_list[[i]], features = c("rna.size", "percent.mt"), ncol = 3) 
  print(qc_plot)
  ggsave(paste0(home_dir, "/Human_Viz/", i, "_qc.pdf"))
} 


#10% for mitochondrial filtering threshold
#+/- 3 median absolute deviations from the median rna library size



rna_size_min <- list() 
rna_size_max <- list() 


for(i in names){ 
  rna_size_min[[i]] <- median(data_list[[i]]@meta.data$rna.size) - (3*mad(data_list[[i]]@meta.data$rna.size)) 
  
  rna_size_max[[i]] <- median(data_list[[i]]@meta.data$rna.size) + (3*mad(data_list[[i]]@meta.data$rna.size)) 
} 



qc_metrics <- c() 
for(i in names){ 
  preqc_cells <- ncol(data_list[[i]]) 
  qc_metrics <- rbind(qc_metrics, preqc_cells) 
}

colnames(qc_metrics) <- "preqc_cells" 
row.names(qc_metrics) <- NULL 
preqc_metrics <- as.data.frame(qc_metrics)



for(i in names){ 
  temp_max <- rna_size_max[[i]] 
  temp_min <-  rna_size_min[[i]] 
  data_list[[i]] <- subset(data_list[[i]], subset = percent.mt < 10 & rna.size < temp_max & rna.size > temp_min) 
  
} 



qc_metrics <- c() 
for(i in names){ 
  preqc_cells <- ncol(data_list[[i]]) 
  qc_metrics <- rbind(qc_metrics, preqc_cells) 
} 



colnames(qc_metrics) <- "postqc_cells" 
row.names(qc_metrics) <- NULL 
postqc_metrics <- as.data.frame(qc_metrics)
cellcounts <- cbind(preqc_metrics, postqc_metrics)
write.csv(cellcounts, file = paste0(home_dir, "/human_cellcounts.csv"))


# Cell Cycle Scoring Addition

s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes 

for(i in names){ 
  data_list[[i]] <- NormalizeData(data_list[[i]]) 
  data_list[[i]] <- CellCycleScoring(data_list[[i]], s.features = s.genes, g2m.features = g2m.genes) 
} 

data_list$AD1@meta.data$Genotype <- "AD"
data_list$AD2@meta.data$Genotype <- "AD"
data_list$C1@meta.data$Genotype <- "C"
data_list$C2@meta.data$Genotype <- "C"


# SCTransform

#Applies SCTransform individually to each sample 

data_list <- lapply(X = data_list, FUN = SCTransform, method = "glmGamPoi", return.only.var.genes = FALSE, vars.to.regress = c("S.Score", "G2M.Score")) 






saveRDS(data_list, file = paste0(home_dir, "/AD_human_objects/combined.RData"))
#data_list <- readRDS(file = paste0(home_dir, "/AD_human_objects/combined.RData"))



variable.features <- SelectIntegrationFeatures(object.list = data_list, nfeatures = 3000) 



# Batch

data_list <- merge(x = data_list[[1]], y = data_list[2:length(data_list)], merge.data=TRUE) 

VariableFeatures(data_list) <- variable.features



data_list <- RunPCA(data_list, verbose = FALSE) 



#Checking differences in principal components 

options(repr.plot.height = 12, repr.plot.width = 20) 
p1 <- DimPlot(object = data_list, reduction = "pca", pt.size = .1, group.by = "orig.ident") 
p2 <- VlnPlot(object = data_list, features = "PC_1", group.by = "orig.ident", pt.size = 0) 
plot_grid(p1,p2) 
ggsave(paste0(home_dir, "/Human_Viz/pca_prebatch.pdf"))

#Batch correction necessary in this instance

#Running Harmony  
options(repr.plot.height = 8, repr.plot.width = 16) 

data_list <- RunHarmony(data_list, assay.use = "SCT", group.by.vars = "orig.ident", plot_convergence = TRUE) 
ggsave(paste0(home_dir, "/Human_Viz/convergence_plot.pdf"))


#Visualizing batch corrected data 

p1 <- DimPlot(object = data_list, reduction = "harmony", pt.size = .1, group.by = "orig.ident") 

p2 <- VlnPlot(object = data_list, features = "harmony_1", group.by = "orig.ident", pt.size = 0) 

plot_grid(p1,p2) 
ggsave(paste0(home_dir, "/Human_Viz/pca_postbatch.pdf"))


#Looking better.

ElbowPlot(data_list, reduction = "harmony", n = 50) 
ggsave(paste0(home_dir, "/Human_Viz/elbowplot.pdf"))


#Lets try 30 PCs

# Reducing Dimensionality
#Performing Dimensional Reduction 

data_list <- RunUMAP(data_list, reduction = "harmony", dims = 1:30, verbose = FALSE) 

data_list <- FindNeighbors(data_list, reduction = "harmony", dims = 1:30, verbose = FALSE) 

data_list <- FindClusters(data_list, resolution = 0.6, verbose = FALSE) 



dittoDimPlot(data_list, "seurat_clusters", do.label = TRUE, legend.show = FALSE, labels.size = 3) 
ggsave(paste0(home_dir, "/Human_Viz/seurat_clusters.pdf"))

dittoDimPlot(data_list, "Genotype") 
ggsave(paste0(home_dir, "/Human_Viz/genotype.pdf"))

dittoDimPlot(data_list, "orig.ident") 
ggsave(paste0(home_dir, "/Human_Viz/orig_ident.pdf"))

dittoDimPlot(data_list, "percent.mt") 
ggsave(paste0(home_dir, "/Human_Viz/percent_mt.pdf"))

dittoDimPlot(data_list, "Phase") 
ggsave(paste0(home_dir, "/Human_Viz/phase.pdf"))




saveRDS(data_list, file = paste0(home_dir, "/AD_human_objects/human_adc_final.RData"))
#data_list <- readRDS(file = paste0(home_dir, "/AD_human_objects/human_adc_final.RData"))

marker_list <- c("UGT8", "KLK6", "KCNH8", "ERMN", "OPALIN",
                 "SLC14A1", "GLIS3", "GLI3", "CHRDL1", "NWD1",
                 "GPR183", "CCL4", "LAPTM5", "CSF1R", "CD14",
                 "APOLD1", "TM4SF1", "FLT1", "A2M",
                 "PDGFRA", "LHFPL3", "MEGF11", "PCDH15",
                 "RELN", "TMEM130", "DLX5", "SST", "EGFR")

for(i in marker_list){
  dittoDimPlot(data_list, var = i, assay = "RNA")
  ggsave(paste0(home_dir, "/Human_Viz/marker_genes/",i, "_umap.pdf"))
  print(VlnPlot(data_list, i, assay = "RNA"))
  ggsave(paste0(home_dir, "/Human_Viz/marker_genes/",i, "_vln.pdf"))
}


# Finding Cluster Identifications
#markers <- FindAllMarkers(data_list, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.5, test.use = "MAST", latent.vars = "orig.ident") 

#clust_counts <- table(data_list@meta.data[,"seurat_clusters"]) 

#clust_markers_list <- list() 

#for(i in unique(markers$cluster)){ 
  
  #clustname <- paste0('cluster', i, '_N', clust_counts[i]) 
  
  #clust_markers_list[[clustname]] <- subset(markers, markers$cluster == i) 
  
  #rownames(clust_markers_list[[clustname]]) <- clust_markers_list[[clustname]]$gene 
  
#} 

#write.xlsx(clust_markers_list, file = paste0(home_dir, "/all_cell_markers_human.xlsx"))

