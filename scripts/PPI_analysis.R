rm(list = ls())

library(Seurat)
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(igraph)

## Declare functions
get_proportion_plot <- function(modules,vec,cols = NULL)
{
  prmat = overlap_modules_wt_categ(mods = modules,vec = vec)
  prdf = melt(prmat);colnames(prdf) = c("X1","X2","value")
  prdf$X1 = factor(prdf$X1)
  prdf$X2 = factor(prdf$X2,levels = names(cols))
  pobj = ggplot(data = prdf,aes(x = X1,y = value,fill = X2)) + geom_bar(stat = "identity",position = "stack") 
  if (!is.null(cols)) {pobj = pobj + scale_fill_manual(values = cols);prdf$X2 = factor(prdf$X2,levels = names(cols))}
  pobj = pobj + theme_bw() + 
    guides(fill = guide_legend(title = "",ncol = 2)) + 
    theme(axis.title = element_blank(),axis.text = element_text(size = 15),legend.position = "bottom")
  return(pobj)
}
overlap_modules_wt_categ <- function(mods,vec)
{
  tpe = setdiff(unique(vec),NA)
  mat = do.call('rbind',lapply(mods,function(x,y,t) {out = table(y[names(y) %in% x])/length(x);out[match(tpe,names(out))]},y = vec,t = tpe))
  mat[is.na(mat)] = 0;
  #rownames(mat) = gsub("^c1_","M",rownames(mat));
  colnames(mat) = tpe
  
  if (any(rowSums(mat,na.rm = TRUE) < 1))
  {
    mat = cbind(mat,1-rowSums(mat,na.rm = TRUE))
    colnames(mat)[ncol(mat)] = "unassigned"
  }
  return(mat)
}

#####################
root.dir <- getwd();
root.dir <- "/media/won-min/My Passport1/SingleCell/Cancer/Melanoma/Data_for_ZNF180_Manuscript";
setwd(root.dir)

# create result directory
dir.create("Results")
out.dir <- "Results/PPI_Analysis";dir.create(out.dir)

# inpuit files
dapfile = "ZNF180_Sequencing_Data/Consensus_DAP.txt" # Consensus DAP file
degfile = "ZNF180_Sequencing_Data/Consensus_DEG.DESeq2.txt" # Consensus DEG file
encodefile = "ZNF180_peaks.anno.txt" # ENCODE ZNF180 TF binding data
ppifile = "9606.protein.links.detailed.v10.5_annotated.RDS"

### load DAP and DEG signatures
dap.comb = read.delim(file = dapfile,sep = "\t",header = T,stringsAsFactors = F)
deg.comb = read.delim(file = degfile,sep = "\t",header = T,stringsAsFactors = F)

sig.list = list("Peak Loss" = unique(subset(dap.comb,DAP.comb == "DN" & grepl("Promoter",A375__annotation))$A375__SYMBOL),
                "DEG Down" = subset(deg.comb,DEG.comb == "DN")$A375__gene.symbol,
                "Peak Gain" = unique(subset(dap.comb,DAP.comb == "UP" & grepl("Promoter",A375__annotation))$A375__SYMBOL),
                "DEG Up" = subset(deg.comb,DEG.comb == "UP")$A375__gene.symbol)

# coerce into list: intersect Peak loss with down-regulated DEGs, and Peak Gain with up-regulated DEGs
int.lst = list(DN = intersect(subset(dap.comb,DAP.comb == "DN" & grepl("Promoter",A375__annotation))$A375__SYMBOL,subset(deg.comb,DEG.comb == "DN")$A375__gene.symbol),
               UP = intersect(subset(dap.comb,DAP.comb == "UP" & grepl("Promoter",A375__annotation))$A375__SYMBOL,subset(deg.comb,DEG.comb == "UP")$A375__gene.symbol))


### load ENCODE ZNF180 TF binding data
encode.peaks <- read.delim(file = encodefile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
encode.peaks$broad.annotation = gsub("Exon \\((.*)$","Exon",gsub("Intron \\((.*)$","Intron",encode.peaks$annotation))
sig.promoter.genes = subset(encode.peaks,grepl("Promoter",broad.annotation))

### Load PPI file from STRING
ppi.data = readRDS(file = ppifile)
ppi.data = subset(ppi.data,combined_score > 900) # subset for high confidence PPIs

# create node pool, based on consensus DEG
bg = Reduce("union",sig.list[c("DEG Down","DEG Up")])
bg = bg[!is.na(bg)]
deg.vec = deg.comb$DEG.comb[match(bg,deg.comb$A375__gene.symbol)] # annotate DEG direction
dap.vec = dap.comb$DAP.comb[match(bg,dap.comb$SKMEL147__SYMBOL)] # annotate consensus DAP direction

node.dat = data.frame(gene.name = bg,DEG = deg.vec,DAP = dap.vec,ZNF180.TFBS = bg %in% sig.promoter.genes$SYMBOL)

# create PPI network
ppi.regulome = subset(ppi.data,symbol1 %in% bg & symbol2 %in% bg)
g.regulome = graph_from_data_frame(d = ppi.regulome[,c("symbol1","symbol2","combined_score")],directed = F,vertices = node.dat)
g.regulome = simplify(g.regulome,edge.attr.comb = "mean")

# identify communities
out = cluster_walktrap(graph = g.regulome)
lst = split(out$names,out$membership);lst = lst[sapply(lst,length) > 10]
print(sapply(lst,length))
lst.rename = paste0("C",1:length(lst))
names(lst) = lst.rename

# update communities on network
vec = rep(NA,vcount(g.regulome))
for (i in 1:length(lst)) vec[V(g.regulome)$name %in% lst[[i]]] = lst.rename[i]
V(g.regulome)$community = factor(vec)
V(g.regulome)$degree = igraph::degree(g.regulome)

# run pathway enrichment
source("scripts/R_functions/enrichment_functions.v3.R")
bg = union(dap.comb$A375__SYMBOL,deg.comb$A375__gene.symbol)
gmtfiles = list.files(path = "MSigDB",full.names = T)
names(gmtfiles) = gsub("^(.*)/|\\.symbols.gmt","",gmtfiles)
msigdb.res = data.frame()
for (i in 1:length(gmtfiles))
{
  hallmark = lapply(read.geneSet(gmtfiles[i]),function(x) x[-1])
  res = perform.AllPairs.FET(lst,hallmark,background = bg)
  res$database = names(gmtfiles)[i]
  msigdb.res = rbind.data.frame(msigdb.res,res)
}
msigdb.res$corrected.FET.pvalue = p.adjust(msigdb.res$FET_pvalue,"BH")
msigdb.res = msigdb.res[order(msigdb.res$FET_pvalue),]
write.table(msigdb.res,file = paste0(out.dir,"/PPI_analysis.MSigDB.txt"),sep = "\t",row.names = F,col.names = T,quote = F)

# make heatmap
library(dplyr)
top.terms = group_by(msigdb.res,set1_Name) %>% slice_min(order_by = FET_pvalue,n = 5)
pdata = subset(msigdb.res,set2_Name %in% top.terms$set2_Name)

pdata$set2_Name = gsub("^GO_|^REACTOME_|^KEGG_","",pdata$set2_Name)

require(reshape2)
m = acast(data = pdata,formula = set1_Name ~ set2_Name,value.var = "FET_pvalue",fun.aggregate = function(x) min(x,na.rm = T))
m[is.infinite(m)] = 1
m = -log10(m)

hobj = hclust(as.dist(sqrt(2*(1-cor(m)))),"complete")
tobj = hclust(dist(m),"complete")

# annotate subclusters by functional groups
library(ggplot2)
msigdb.heat = ggplot() + 
  geom_tile(data = pdata,
            aes(y = set2_Name,x = set1_Name,fill = -log10(corrected.FET.pvalue))) + 
  scale_fill_gradient2(low = "white",mid = "white",high = "red") + 
  scale_y_discrete(limits = hobj$labels[hobj$order]) + 
  geom_point(data = subset(pdata,corrected.FET.pvalue < 0.05),
             aes(y = set2_Name,x = set1_Name,size = enrichment.foldchange)) + 
  scale_x_discrete(limits = tobj$labels[tobj$order]) + 
  guides(fill = guide_colorbar(title = "-log10(FET FDR)"),size = guide_legend(title = "EFC")) + 
  theme_classic() + 
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 14,angle = 45,vjust= 1,hjust =1),
        axis.text.y = element_text(size = 14),legend.title = element_text(size = 15),legend.text = element_text(size = 13),
        strip.text.x.top = element_text(angle = 90,hjust = 0.5,size = 14),
        legend.position = "bottom",legend.direction = "horizontal")

#png(file = paste0(out.dir,"/PPI_analysis.MSigDB.png"),res = 450,width = 6500,height = 5500)
print(msigdb.heat)
#dev.off()

### name communities by functions
function.name = c("C7" = "RNA binding protein/Ubiquitins","C2" = "Protein folding","C3" = "NRAS/MAPK/RNA Splicing","C1" = "FA/Choromosome organizaiton",
                  "C4" = "DNA Repair/Replication","C5" = "Endocytosis","C6" = "Chaperone","C8" = "Spliceosome/HS-GAG biosynthesis","C9" = "MHC-I Antigen presentation/Proteasome")

# check ratio of DEGs
node.dat$community = V(g.regulome)$community[match(node.dat$gene.name,V(g.regulome)$name)]
node.dat$community.function = paste0(node.dat$community,":",function.name[node.dat$community])

pdat = subset(node.dat,!is.na(community))
vec = pdat$DEG;names(vec) = pdat[[1]]
out = get_proportion_plot(modules = split(pdat$gene.name,factor(pdat$community.function)),vec = vec,cols = c("UP" = "red","DN" = "blue")) +
  scale_y_continuous(labels = scales::percent) + 
  guides(fill = "none") + theme_minimal() + 
  theme(axis.title = element_blank(),axis.text.y = element_text(size = 15),axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1)) +
  coord_flip()

print(out)