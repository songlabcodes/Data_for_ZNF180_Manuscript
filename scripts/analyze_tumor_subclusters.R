rm(list = ls())

library(ggplot2)
library(cowplot)
library(Seurat)

### run DEG per cell cluster
library(DESeq2)
library(Matrix.utils)
library(Matrix)
library(doParallel)

#######
root.dir <- getwd();
setwd(root.dir)
out.dir <- "Results/Jerby_Arnon_et_al_2018"

dir.create("Results")
dir.create(out.dir)

tumorfile = "Jerby_Arnon_et_al_2018/integrated.tumor.RDS"

### load single cell data
seu.tumor = readRDS(file = tumorfile)

#### Find immune checkpoint genes in tumors with ZNF180 pathways
if (TRUE)
{
  add_axis <- function(pobj,xy)
  {
    out = pobj + annotate(geom = "segment",x = min(xy[,1]),xend = (max(xy[,1])-min(xy[,1]))*0.25 + min(xy[,1]),y = min(xy[,2]),yend = min(xy[,2]),colour = "black",arrow = grid::arrow(length = unit(0.015, "npc"))) + 
      annotate(geom = "segment",x = min(xy[,1]),xend = min(xy[,1]),y = min(xy[,2]),yend = (max(xy[,2])-min(xy[,2]))*0.25 + min(xy[,2]),colour = "black",arrow = grid::arrow(length = unit(0.015, "npc")) )+
      annotate(geom = "text",x = mean(c(min(xy[,1]),(max(xy[,1])-min(xy[,1]))*0.25 + min(xy[,1]))),y = min(xy[,2]),hjust = 0.4,vjust = -0.5,label = "UMAP1") + 
      annotate(geom = "text",x = min(xy[,1]),y = mean(c(min(xy[,2]),(max(xy[,2])-min(xy[,2]))*0.25 + min(xy[,2]))),angle = 90,hjust = 0.5,vjust = 1.5,label = "UMAP2")
    return(out)
  }
  
  gene.name = rownames(seu.tumor)
  genes = gene.name[grepl("^AXL$|PDCD1LG2|CD274|HLA-A|HLA-B|HLA-C|^PVR$|PVRL[0-9]$|NECTIN|NECTIN2|ZNF180|^MSH2$|^FANCM$|^ATR$|^ATM$|MITF|HLA-DR|HLA-DP|HLA-DQ|IFNG",gene.name)]
  
  # mark gene category
  vec = rep(NA,length(genes))
  vec[grepl("HLA-A|HLA-B|HLA-C",genes)] = "MHC-I"
  vec[grepl("HLA-DR|HLA-DP|HLA-DQ",genes)] = "MHC-II"
  vec[grepl("CD274|PDCD1LG2",genes)] = "PD-L1/-L2"
  vec[grepl("ZNF180|MSH2|FANCM|ATR|ATM",genes)] = "ZNF180\ndriven\nDDR"
  vec[grepl("MITF|AXL",genes)] = "MITF/\nAXL"
  vec[grepl("PVR|PVRL|NECTIN|PVRL",genes)] = "CD155/\nNECTIN"
  vec[grepl("IFNG",genes)] = "IFN-\u03B3"
  
  gene.annot = data.frame(genes = genes,category = vec)
  
  gene.vec = genes;names(gene.vec) = vec
  dot.obj = DotPlot(seu.tumor,features = gene.vec,assay = "SCT") + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
    theme(legend.position = "bottom",legend.direction = "horizontal",strip.text.x.top = element_text(angle = 90,hjust = 0.5,size = 15),
          axis.title = element_blank())
  
  xy= Embeddings(seu.tumor,"umap")
  umap.obj.cls = DimPlot(seu.tumor,label = T,label.size = 8) + guides(colour = "none") + 
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),
          plot.background = element_blank(),plot.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  umap.obj.cls = add_axis(umap.obj.cls,xy)
  
  umap.obj.trt = DimPlot(seu.tumor,group.by = "treatment.group") + 
    scale_colour_manual(values = c("post.treatment" = "brown","treament.naive" = "chartreuse")) + 
    guides(colour = "none") + 
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),
          plot.background = element_blank(),plot.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  umap.obj.trt = add_axis(umap.obj.trt,xy)
  
  # show expressions of key genes
  library(cowplot)
  plot_grid(plot_grid(umap.obj.cls,umap.obj.trt,ncol = 2,labels = c("A","B"),label_size = 30),dot.obj,ncol = 1,rel_heights = c(0.4,0.6),labels = c("","C"),label_size = 30)
  
}
