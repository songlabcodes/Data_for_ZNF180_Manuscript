rm(list = ls())

library(Seurat)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(UCell)
library(ggvenn)
library(cowplot)

root.dir <- getwd()
setwd(root.dir)
seuratfile = "Wouters_et_al_2020/integrated.seurat.RDS" # seurat object
dapfile = "ZNF180_Sequencing_Data/Consensus_DAP.txt" # Consensus DAP file
degfile = "ZNF180_Sequencing_Data/Consensus_DEG.DESeq2.txt" # Consensus DEG file

# load the single cell transcriptome
seu = readRDS(file = seuratfile)
seu@meta.data$cell.line = gsub("^(.*)10x_|_BL$","",seu@meta.data$data.id)

# load ATAC and RNA seq data
dap.comb = read.delim(file = dapfile,sep = "\t",header = T,stringsAsFactors = F)
deg.comb = read.delim(file = degfile,sep = "\t",header = T,stringsAsFactors = F)

# create Venn diagram

sig.list = list("Peak Loss" = unique(subset(dap.comb,DAP.comb == "DN" & grepl("Promoter",A375__annotation))$A375__SYMBOL),
                "DEG Down" = subset(deg.comb,DEG.comb == "DN")$A375__gene.symbol,
                "Peak Gain" = unique(subset(dap.comb,DAP.comb == "UP" & grepl("Promoter",A375__annotation))$A375__SYMBOL),
                "DEG Up" = subset(deg.comb,DEG.comb == "UP")$A375__gene.symbol)

print(ggvenn(sig.list,fill_color = c("deepskyblue","blue","coral1","red")))

# coerce into list: List of epigenetically regulated genes by ZNF180
int.lst = list(DN = intersect(subset(dap.comb,DAP.comb == "DN" & grepl("Promoter",A375__annotation))$A375__SYMBOL,subset(deg.comb,DEG.comb == "DN")$A375__gene.symbol),
               UP = intersect(subset(dap.comb,DAP.comb == "UP" & grepl("Promoter",A375__annotation))$A375__SYMBOL,subset(deg.comb,DEG.comb == "UP")$A375__gene.symbol))
names(int.lst) = paste0("ZNF180.",names(int.lst))

# get module score per cell
seu <- AddModuleScore_UCell(seu,features=int.lst, name="_UCell")

cell.order = c("MM001","MM011","MM031","MM074","MM087","MM057","A375","MM029","MM047","MM099")
cell.state = list("Melanocytic" = c("MM001","MM011","MM031"),
                  "Intermediate" = c("MM074","MM087","MM057"),
                  "Mesenchymal-like" = c("MM029","MM047","MM099"))

vec = rep(ncol(seu))
for (i in 1:length(cell.state))
{
  vec[seu@meta.data$cell.line %in% cell.state[[i]]] = names(cell.state)[i]
}
seu@meta.data$cell.state = vec;

#### make plot
# get data
a = DotPlot(seu,features = c("ZNF180.DN_UCell","ZNF180.UP_UCell"),
            group.by = "cell.line",cluster.idents = F,scale = T) + 
  scale_colour_gradient2(low = "blue",mid = "white",high = "red",midpoint = 0) +
  guides(size = "none",colour = guide_colorbar(title = "Scaled\nScore",title.position = "top")) + 
  scale_size(range = c(6,6)) +
  labs(title = "shZNF180\nUCell Score") + 
  scale_y_discrete(limits = cell.order) + 
  scale_x_discrete(labels = c("DN","UP")) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),axis.title = element_blank(),legend.position = "bottom",legend.direction = "vertical")

b = DotPlot(seu,features = c("MITF","AXL","ZNF180","FOSL1","HMGA1","JUN","PVR"),
            group.by = "cell.line") + scale_size(range = c(2,6)) +
  scale_y_discrete(limits = cell.order) +
  labs(title = "Cell State Markers") +
  guides(colour = guide_colorbar(title = "Avg.\nExpression",title.position = "top",legend.direction = "vertical"),size = guide_legend(title = "%.\nExpressed",nrow = 2,title.position = "top")) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),axis.title = element_blank(),legend.position = "bottom",legend.direction = "vertical",
        #axis.text.y = element_blank(),
        plot.title = element_text(hjust= 0.5))

plot_grid(a,b,ncol = 2,rel_widths = c(0.3,0.7),axis = "tblr",align = "h")
