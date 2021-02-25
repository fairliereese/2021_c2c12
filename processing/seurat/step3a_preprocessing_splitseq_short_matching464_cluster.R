# Process 1k lib Seurat object and make plots

library(dplyr)
library(Seurat)
library(Matrix)
library(stringr)
library(Matrix.utils)
library(ggplot2, quietly=T)

setwd("/share/crsp/lab/seyedam/share/c2c12_pacbio_singlecell/scRNA_illumina/")

load("seurat_objects/mb_mt_464_raw_seuratobj.rda")

#################### Process 464 cells ####################
mb_mt_1k <- SCTransform(mb_mt_1k, 
                        vars.to.regress = c("percent.mt","nFeature_RNA"), verbose = TRUE)
mb_mt_1k <- RunPCA(mb_mt_1k, features = VariableFeatures(object = mb_mt_1k)) # 18,313
mb_mt_1k <- FindNeighbors(mb_mt_1k, dims = 1:15)
mb_mt_1k <- FindClusters(mb_mt_1k, resolution = 0.9,algorithm = 4)
mb_mt_1k <- RunUMAP(mb_mt_1k, dims = 1:15)

mb_mt_1k@meta.data$final_clusters_ordered  = mb_mt_1k@meta.data$seurat_clusters 
Idents(mb_mt_1k) = mb_mt_1k@meta.data$final_clusters_ordered

# Rename clusters to be in rough order of differentiation
mb_mt_1k=RenameIdents(mb_mt_1k,
                      '1'='2', 
                      '2'='5',
                      '3'='1',
                      '4'='7',
                      '5'='3',
                      '6'='4',
                      '7'='6')
mb_mt_1k@meta.data$final_clusters_ordered = Idents(mb_mt_1k)
my_levels <- c("1","2","3","4","5","6","7")
mb_mt_1k@meta.data$final_clusters_ordered <- factor(mb_mt_1k@meta.data$final_clusters_ordered, levels= rev(my_levels))

cellID_clusters  = as.data.frame(mb_mt_1k$final_clusters_ordered)
write.table(cellID_clusters,file="shortreads_464cells_7clustersOrdered_IDs.tsv", quote = F)
save(mb_mt_1k,file="seurat_objects/mb_mt_1k_processed.rda")


########
########