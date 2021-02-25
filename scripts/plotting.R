library(Seurat)
library(tidyverse)
options(stringsAsFactors = FALSE)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

get_sample_colors <- function() {
    
    mt_green = '#019f73'
    mb_pink = '#cb79a7'
    mb_blue = '#57b4e9'
    colors = list(mb_pink, mb_blue, mt_green)
    names(colors) = c('MB_cells', 'MB_nuclei', 'MT_nuclei')
    
    return(colors)
}

get_long_clust_colors <- function() {
    purple = '#8DA0CB'
    yellow = '#FFD92F'
    green = '#A6D854'
    pink = '#E78AC3'
    beige = '#E5C494'
    salmon = '#FC8D62'
    teal = '#66C2A5'
    mint = '#b4fade'

    long_clusters = list(1, 2, 3, 4, 5, 6, 7)
    colors = c(teal, purple, yellow, green, beige, salmon, pink)
    
    return(colors)
}

plot_clust_prop <- function(seurat_obj, gb, colors, opref) {
    clusters = unique(seurat_obj$leiden)
    clusters = clusters[order(clusters)]
    
    prop_table = data.frame()
    for(clust in clusters){
        curr = subset(seurat_obj@meta.data, leiden == clust)
        temp = as.data.frame(table(curr[,gb])/nrow(curr))
        temp$leiden = clust
        prop_table = rbind(prop_table, temp)
    }
    prop_table$leiden = as.character(prop_table$leiden)
    prop_table$leiden <- factor(prop_table$leiden, levels=c(7,6,5,4,3,2,1))
    
    p = ggplot(prop_table, aes(y=Freq, x=leiden, fill=Var1))+
    geom_bar(stat='identity')+
    scale_fill_manual(values=colors)+
    coord_flip()+
    ylab('Fraction')
    
    fname = paste('figures/gene', gb, 'cluster_props.pdf', sep='_')
    pdf(fname, width=.75*7, height=2.75*5)
    print(p)
    dev.off()
    
    return(p)
}

plot_short_clust_prop <- function(seurat_obj, gb, colors, opref) {
    clusters = unique(seurat_obj$short_leiden)
    clusters = clusters[order(clusters)]
    
    prop_table = data.frame()
    for(clust in clusters){
        curr = subset(seurat_obj@meta.data, short_leiden == clust)
        temp = as.data.frame(table(curr[,gb])/nrow(curr))
        temp$short_leiden = clust
        prop_table = rbind(prop_table, temp)
    }
    prop_table$short_leiden = as.character(prop_table$short_leiden)
    prop_table$short_leiden <- factor(prop_table$short_leiden, levels=c(7,6,5,4,3,2,1))
    
    p = ggplot(prop_table, aes(y=Freq, x=short_leiden, fill=Var1))+
    geom_bar(stat='identity')+
    scale_fill_manual(values=colors)+
    coord_flip()+
    ylab('Fraction')
    
    fname = paste('figures/gene', gb, 'short_cluster_props.pdf', sep='_')
    pdf(fname, width=.75*7, height=2.75*5)
    print(p)
    dev.off()
    
    return(p)
}