# Make Seurat object, QC filter for Velocyto

library(dplyr)
library(Seurat)
library(Matrix)
library(stringr)
library(Matrix.utils)
library(ggplot2, quietly=T)

setwd("/share/crsp/lab/seyedam/share/c2c12_pacbio_singlecell/scRNA_illumina/")

#################### Make separate Seurat objects for all libs ####################
load("splitbio_pipeline_output/mb_mt_counts_all.rda")
load("splitbio_pipeline_output/metadata_all.rda")

# Seurat takes matrix of genes in rows and cells in columns 
mb_mt <- CreateSeuratObject(counts = mb_mt_counts_all, project = "C2C12", min.cells = 0, min.features = 0)
mb_mt@meta.data$Library <- sapply(strsplit(as.character(colnames(mb_mt)), "_"), "[[", 2)
mb_mt@meta.data$SampleType <- metadata_all$sample
mb_mt@meta.data$CellType <-  sapply(strsplit(as.character(metadata_all$sample), "_"), "[[", 1)
table(mb_mt@meta.data$CellType) # 18,406 MB and 20,987 MT
table(mb_mt@meta.data$SampleType) # 7,960 MB cells and 10,446 MB nuclei
mb_mt[["percent.mt"]] <- PercentageFeatureSet(mb_mt, pattern = "mt-")

#################### Seurat QC cutoffs for big object ####################
mb_mt_filt <- subset(
  x = mb_mt,
  subset = 
    nFeature_RNA > 500 &
    nCount_RNA < 200000 &
    percent.mt < 20
)

table(mb_mt_filt$SampleType) # 7,797 MB cells, 10,194 MB nuc, 18,878 MT nuc
save(mb_mt_filt, file="seurat_objects/mb_mt_filt_raw_seuratobj.rda")

#################### Make separate Seurat objects for 1k lib ####################
load("splitbio_pipeline_output/mb_mt_counts_1k_lib.rda")
load("splitbio_pipeline_output/metadata_1k.rda")

# Make sure all of the matching long read cells are in the final data
longread_bcs = read.table("ref/bcs_for_liz.txt", stringsAsFactors = F)
table(longread_bcs$V1 %in% colnames(mb_mt_filt)) # yes

# Use 464 cells that passed long read filters
mb_mt_counts_1k_464 = mb_mt_counts_1k[,colnames(mb_mt_counts_1k) %in% paste0(longread_bcs$V1)]
metadata_1k_464 = metadata_1k[metadata_1k$cell_barcode_24nt %in% longread_bcs$V1,]
metadata_1k_464 = metadata_1k_464[match(colnames(mb_mt_counts_1k_464), metadata_1k_464$cell_barcode_24nt),]

table(metadata_1k_464$cell_barcode_24nt == colnames(mb_mt_counts_1k_464))

mb_mt_1k <- CreateSeuratObject(counts = mb_mt_counts_1k_464, project = "C2C12", min.cells = 0, min.features = 0)
mb_mt_1k@meta.data$Library <- sapply(strsplit(as.character(colnames(mb_mt_1k)), "_"), "[[", 2)
mb_mt_1k@meta.data$SampleType <- metadata_1k_464$sample
mb_mt_1k@meta.data$CellType <-  sapply(strsplit(as.character(metadata_1k_464$sample), "_"), "[[", 1)
table(mb_mt_1k@meta.data$CellType) # 255 MB 209 MT
table(mb_mt_1k@meta.data$SampleType) # 110 MB cells and 145 MB nuclei
mb_mt_1k[["percent.mt"]] <- PercentageFeatureSet(mb_mt_1k, pattern = "mt-")

save(mb_mt_1k, file="seurat_objects/mb_mt_464_raw_seuratobj.rda")

########
# write barcodes for velocyto 
write.csv(colnames(mb_mt_filt[["RNA"]]),file="ref/velocyto_barcodes_expt2_9kand1k_final.csv",quote=F, row.names = F)

metadata_all_filt = metadata_all[metadata_all$cell_barcode_24nt %in% colnames(mb_mt_filt[["RNA"]]),]

metadata_final = metadata_all_filt[match(colnames(mb_mt_filt[["RNA"]]), metadata_all_filt$cell_barcode_24nt),]
table(metadata_final$cell_barcode_24nt %in% colnames(mb_mt_filt[["RNA"]]))

write.csv(metadata_final,file="ref/metadata_final_36869cells.csv", quote = F, row.names = F)

########
########