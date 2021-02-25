# Read in Splitseq data to make Seurat object

library(plyr)
library(dplyr)
library(Matrix)
library(stringr)
library(Matrix.utils)
library(Seurat)

setwd("/share/crsp/lab/seyedam/share/c2c12_pacbio_singlecell/scRNA_illumina/")

options("stringsAsFactors" = FALSE)

#################### Read in barcode 1 sequences #################### 
rnd1_barcode = as.data.frame(read.csv('ref/bc_8nt_v2.csv'))

#################### 9kA #################### 
celldata_9kA = read.csv('splitbio_pipeline_output/9kA/cell_metadata.csv')
counts_9kA = readMM('splitbio_pipeline_output/9kA/DGE.mtx')
genelist_9kA = read.csv('splitbio_pipeline_output/9kA/genes.csv')
colnames(counts_9kA) = paste0(genelist_9kA$gene_id,"_",genelist_9kA$gene_name)

# get 18nt sequence from rounds 3 and 2
celldata_9kA$rnd3_rnd2 = do.call("rbind", strsplit(as.character(celldata_9kA$cell_barcode), "_"))[,1]
celldata_9kA$Index = do.call("rbind", strsplit(as.character(celldata_9kA$cell_barcode), "_"))[,2]

# get 6nt sequence from round 1 to make 24nt barcode 
celldata_9kA2 = join(celldata_9kA,rnd1_barcode)
celldata_9kA2$cell_barcode_24nt = paste0(celldata_9kA2$rnd3_rnd2,celldata_9kA2$Barcode,"_9kA")
rownames(counts_9kA) = celldata_9kA2$cell_barcode_24nt

#################### 9kB #################### 
celldata_9kB = read.csv('splitbio_pipeline_output/9kB/cell_metadata.csv')
counts_9kB = readMM('splitbio_pipeline_output/9kB/DGE.mtx')
genelist_9kB = read.csv('splitbio_pipeline_output/9kB/genes.csv')
colnames(counts_9kB) = paste0(genelist_9kB$gene_id,"_",genelist_9kB$gene_name)

# get 18nt sequence from rounds 3 and 2
celldata_9kB$rnd3_rnd2 = do.call("rbind", strsplit(as.character(celldata_9kB$cell_barcode), "_"))[,1]
celldata_9kB$Index = do.call("rbind", strsplit(as.character(celldata_9kB$cell_barcode), "_"))[,2]

# get 6nt sequence from round 1 to make 24nt barcode 
celldata_9kB2 = join(celldata_9kB,rnd1_barcode)
celldata_9kB2$cell_barcode_24nt = paste0(celldata_9kB2$rnd3_rnd2,celldata_9kB2$Barcode,"_9kB")
rownames(counts_9kB) = celldata_9kB2$cell_barcode_24nt

#################### 9kC #################### 
celldata_9kC = read.csv('splitbio_pipeline_output/9kC/cell_metadata.csv')
counts_9kC = readMM('splitbio_pipeline_output/9kC/DGE.mtx')
genelist_9kC = read.csv('splitbio_pipeline_output/9kC/genes.csv')
colnames(counts_9kC) = paste0(genelist_9kC$gene_id,"_",genelist_9kC$gene_name)

# get 18nt sequence from rounds 3 and 2
celldata_9kC$rnd3_rnd2 = do.call("rbind", strsplit(as.character(celldata_9kC$cell_barcode), "_"))[,1]
celldata_9kC$Index = do.call("rbind", strsplit(as.character(celldata_9kC$cell_barcode), "_"))[,2]

# get 6nt sequence from round 1 to make 24nt barcode 
celldata_9kC2 = join(celldata_9kC,rnd1_barcode)
celldata_9kC2$cell_barcode_24nt = paste0(celldata_9kC2$rnd3_rnd2,celldata_9kC2$Barcode,"_9kC")
rownames(counts_9kC) = celldata_9kC2$cell_barcode_24nt

#################### 9kD #################### 
celldata_9kD = read.csv('splitbio_pipeline_output/9kD/cell_metadata.csv')
counts_9kD = readMM('splitbio_pipeline_output/9kD/DGE.mtx')
genelist_9kD = read.csv('splitbio_pipeline_output/9kD/genes.csv')
colnames(counts_9kD) = paste0(genelist_9kD$gene_id,"_",genelist_9kD$gene_name)

# get 18nt sequence from rounds 3 and 2
celldata_9kD$rnd3_rnd2 = do.call("rbind", strsplit(as.character(celldata_9kD$cell_barcode), "_"))[,1]
celldata_9kD$Index = do.call("rbind", strsplit(as.character(celldata_9kD$cell_barcode), "_"))[,2]

# get 6nt sequence from round 1 to make 24nt barcode 
celldata_9kD2 = join(celldata_9kD,rnd1_barcode)
celldata_9kD2$cell_barcode_24nt = paste0(celldata_9kD2$rnd3_rnd2,celldata_9kD2$Barcode,"_9kD")
rownames(counts_9kD) = celldata_9kD2$cell_barcode_24nt

#################### 9kE #################### 
celldata_9kE = read.csv('splitbio_pipeline_output/9kE/cell_metadata.csv')
counts_9kE = readMM('splitbio_pipeline_output/9kE/DGE.mtx')
genelist_9kE = read.csv('splitbio_pipeline_output/9kE/genes.csv')
colnames(counts_9kE) = paste0(genelist_9kE$gene_id,"_",genelist_9kE$gene_name)

# get 18nt sequence from rounds 3 and 2
celldata_9kE$rnd3_rnd2 = do.call("rbind", strsplit(as.character(celldata_9kE$cell_barcode), "_"))[,1]
celldata_9kE$Index = do.call("rbind", strsplit(as.character(celldata_9kE$cell_barcode), "_"))[,2]

# get 6nt sequence from round 1 to make 24nt barcode 
celldata_9kE2 = join(celldata_9kE,rnd1_barcode)
celldata_9kE2$cell_barcode_24nt = paste0(celldata_9kE2$rnd3_rnd2,celldata_9kE2$Barcode,"_9kE")
rownames(counts_9kE) = celldata_9kE2$cell_barcode_24nt

#################### 9kF #################### 
celldata_9kF = read.csv('splitbio_pipeline_output/9kF/cell_metadata.csv')
counts_9kF = readMM('splitbio_pipeline_output/9kF/DGE.mtx')
genelist_9kF = read.csv('splitbio_pipeline_output/9kF/genes.csv')
colnames(counts_9kF) = paste0(genelist_9kF$gene_id,"_",genelist_9kF$gene_name)

# get 18nt sequence from rounds 3 and 2
celldata_9kF$rnd3_rnd2 = do.call("rbind", strsplit(as.character(celldata_9kF$cell_barcode), "_"))[,1]
celldata_9kF$Index = do.call("rbind", strsplit(as.character(celldata_9kF$cell_barcode), "_"))[,2]

# get 6nt sequence from round 1 to make 24nt barcode 
celldata_9kF2 = join(celldata_9kF,rnd1_barcode)
celldata_9kF2$cell_barcode_24nt = paste0(celldata_9kF2$rnd3_rnd2,celldata_9kF2$Barcode,"_9kF")
rownames(counts_9kF) = celldata_9kF2$cell_barcode_24nt

#################### 1k #################### 
celldata_1k = read.csv('splitbio_pipeline_output/1k/cell_metadata.csv')
counts_1k = readMM('splitbio_pipeline_output/1k/DGE.mtx')
genelist_1k = read.csv('splitbio_pipeline_output/1k/genes.csv')
colnames(counts_1k) = paste0(genelist_1k$gene_id,"_",genelist_1k$gene_name)

# get 18nt sequence from rounds 3 and 2
celldata_1k$rnd3_rnd2 = do.call("rbind", strsplit(as.character(celldata_1k$cell_barcode), "_"))[,1]
celldata_1k$Index = do.call("rbind", strsplit(as.character(celldata_1k$cell_barcode), "_"))[,2]

# get 6nt sequence from round 1 to make 24nt barcode 
celldata_1k2 = join(celldata_1k,rnd1_barcode)
celldata_1k2$cell_barcode_24nt = paste0(celldata_1k2$rnd3_rnd2,celldata_1k2$Barcode)
rownames(counts_1k) = celldata_1k2$cell_barcode_24nt

#################### Combine sublibrary matrices by matching ensembl ID #################### 
counts = RowMergeSparseMatrices(t(counts_9kA), t(counts_9kB))
dim(counts)

counts = RowMergeSparseMatrices(counts, t(counts_9kC))
dim(counts)

counts = RowMergeSparseMatrices(counts, t(counts_9kD))
dim(counts)

counts = RowMergeSparseMatrices(counts, t(counts_9kE))
dim(counts)

counts = RowMergeSparseMatrices(counts, t(counts_9kF))
dim(counts)

#################### Make counts matrices for all 7 merged libs and 1k lib on its own #################### 
# 1 1k sublibrary
# 20,557 genes, 568 cells
mb_mt_counts_1k = t(as.matrix(counts_1k))
rownames(mb_mt_counts_1k) = sapply(strsplit(as.character(rownames(mb_mt_counts_1k)), "_"), "[[", 2)

# Merge 9k libs and 1k for full matrix
# 24,373 genes, 39,393 cells
mb_mt_counts_all = RowMergeSparseMatrices(counts, t(counts_1k))

#################### Make metadata files for 6 combined 9k sublibraries, 1k by himself, and all 7 merged #################### 
# Merge 9k metadata files
metadata_9k= rbind(celldata_9kA2,celldata_9kB2,celldata_9kC2,celldata_9kD2,celldata_9kE2,celldata_9kF2)

# Check that 24nt cell barcode in metadata (plus _sublibrary ID) matches column names
table(metadata_9k$cell_barcode_24nt == colnames(mb_mt_counts)) 

# Get metadata for 1k library
metadata_1k = celldata_1k2

# Check that 24nt cell barcode in metadata (plus _sublibrary ID) matches column names
table(metadata_1k$cell_barcode_24nt == colnames(mb_mt_counts_1k))

# Merge all metadata
metadata_all= rbind(celldata_9kA2,celldata_9kB2,celldata_9kC2,celldata_9kD2,celldata_9kE2,celldata_9kF2,metadata_1k)
table(metadata_all$cell_barcode_24nt == colnames(mb_mt_counts_all)) 

#################### Save tsv and rda files #################### 
write.table(metadata_1k, file="splitbio_pipeline_output/metadata_1k.tsv", row.names = F, quote = F)
write.table(mb_mt_counts_1k, file="splitbio_pipeline_output/mb_mt_counts_1k_lib.tsv", quote = F)
save(metadata_1k, file="splitbio_pipeline_output/metadata_1k.rda")
save(mb_mt_counts_1k, file="splitbio_pipeline_output/mb_mt_counts_1k_lib.rda")

write.table(metadata_all, file="splitbio_pipeline_output/metadata_all_7_libs.tsv", row.names = F, quote = F)
write.table(mb_mt_counts_all, file="splitbio_pipeline_output/mb_mt_counts_all.tsv", quote = F)
save(metadata_all, file="splitbio_pipeline_output/metadata_all.rda")
save(mb_mt_counts_all, file="splitbio_pipeline_output/mb_mt_counts_all.rda")

#############
#############