library(Seurat)
library(viridis)
library(stringr)

sampletype_cols = c('#cb79a7','#57b4e9','#019f73')
sampletype_cols_mnc = c('#CB79A7','#019F73','#08684C')
cluster20_colors = c("#CBD844","#C193C6","#FAD53E","#E28BC3","#D7BF9C",
                     "#C79693","#66C2A5","#E7C689","#C5B9A7","#9F9BC9",
                     "#AAD852","#9E9DBA","#9DAE8C","#D49A73","#ECD836",
                     "#D2A29F","#BABF77","#B3B3B3","#F08F6D","#F1CD64")

load("seurat_objects/mb_mt_1k_processed.rda")
load("seurat_objects/mb_mt_36869_cells_20clusters.rda")
load("seurat_objects/mb_mt_filt_raw_seuratobj.rda")

#################### Supplemental 2C, violin plots of QC-filtered 464 cells (short reads) ####################
fname = "figures/qc_violinplot_464cells.pdf"
pdf(file = fname,
    width = 6,
    height = 5)
VlnPlot(mb_mt_1k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "SampleType", pt.size = 0.5, ncol = 3, cols = sampletype_cols)
dev.off()

#################### Supplemental 2E, violin plots of all QC-filtered cells (short reads) ####################
# using pre-velocyto filtered object
fname = "figures/qc_violinplot_allcells.pdf"
pdf(file = fname,
    width = 6,
    height = 5)
VlnPlot(mb_mt_filt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "SampleType", pt.size = 0, ncol = 3, cols = sampletype_cols)
dev.off()

#################### Supplemental 2F ####################
genes = c("Mki67","Hmga2","Pax7","Igfbp5","Col3a1","Col1a1","Myog","Mybph","Myh3","Rbm24")
fname = "figures/464_shortreads_featureplots_ordered.pdf"
pdf(file = fname,
    width = 10, 
    height = 4)
FeaturePlot(mb_mt_1k, order=T,ncol=5,features = genes) & NoLegend() & NoAxes() & 
  scale_colour_gradientn(colours = viridis(11))
dev.off()

#################### Supplemental 2G ####################
mb_mt@meta.data$grouped_clusters  = mb_mt@meta.data$final_clusters_ordered 
Idents(mb_mt) = mb_mt@meta.data$grouped_clusters

mb_mt=RenameIdents(mb_mt,
                   '1'='MB', 
                   '2'='MB',
                   '3'='MB',
                   '4'='MB',
                   '5'='MB',
                   '6'='MB',
                   '7'='MB',
                   '8'='MNC',
                   '9'='MNC',
                   '10'='MNC',
                   '11'='MNC',
                   '12'='MNC',
                   '13'='MNC',
                   '14'='MNC',
                   '15'='MNC',
                   '16'='MT',
                   '17'='MT',
                   '18'='MT',
                   '19'='MT',
                   '20'='MT')
mb_mt@meta.data$grouped_clusters = Idents(mb_mt)

fname = "figures/shortread_umap_mb_mt_mnc.pdf"
pdf(file = fname,
    width = 5.5, 
    height = 5)
DimPlot(
  object = mb_mt,
  group.by = "grouped_clusters",
  label = F,label.size = 6,
  repel = TRUE)  + NoLegend() + NoAxes()+ scale_color_manual(values = sampletype_cols_mnc)
dev.off()

#################### Supplemental 2H ####################
cellID_clusters  = as.data.frame(mb_mt$final_clusters_ordered)
cellID_clusters$cellID  =sapply(strsplit( sapply(strsplit(rownames(cellID_clusters), ":"), "[[", 2), "_"), "[[", 1)
cellID_clusters$sublib  =sapply(strsplit( sapply(strsplit(rownames(cellID_clusters), ":"), "[[", 2), "_"), "[[", 2)

# Get cell ID and cluster from big object
cellID_clusters_464 = cellID_clusters[cellID_clusters$sublib == "1kx" & cellID_clusters$cellID %in% colnames(mb_mt_1k),]
cellID_clusters_464 = cellID_clusters_464[match(colnames(mb_mt_1k), cellID_clusters_464$cellID),]
table(cellID_clusters_464$cellID == colnames(mb_mt_1k)) # yes

mb_mt_1k$clusters_20 = cellID_clusters_464$`mb_mt$final_clusters_ordered`

fname = "figures/shortread_umap_464cells_20clust.pdf"
pdf(file = fname,
    width = 5.5, 
    height = 5)
DimPlot(
  object = mb_mt_1k,pt.size = 2.2,
  group.by = "clusters_20",
  label = F,label.size = 6,
  repel = TRUE)  + NoLegend() + NoAxes()+ scale_color_manual(values = cluster20_colors)
dev.off()

#################### Supplemental 2I ####################
Idents(mb_mt) = mb_mt@meta.data$final_clusters_ordered
cluster.averages <- AverageExpression(mb_mt, return.seurat = TRUE)

mb_mt.markers <- FindAllMarkers(mb_mt, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
mb_mt.markers = mb_mt.markers[mb_mt.markers$p_val_adj < 0.01,] # this is the same df as the marker table in the processing directory

#mb_mt.markers = read.table("../processing/tables/seurat/splitseq_markers_minpct0.1_lfc0.1_20Clusters_fdr0.01.tsv", header = T, stringsAsFactors = F)

top10 <- mb_mt.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

fname = "figures/heatmap_20clusters_top10.pdf"
pdf(file=fname,
    width = 8, 
    height = 9.5)
DoHeatmap(cluster.averages, label=T,features = top10$gene, group.colors = cluster20_colors, raster=F,draw.lines = FALSE) + 
  scale_fill_viridis() + theme(axis.text.y = element_text(size = 6))
dev.off()


########
########