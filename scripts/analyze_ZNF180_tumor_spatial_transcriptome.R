rm(list = ls())

#### load necessary libraries and function
library(Giotto)
library(InstaPrism)
library(RColorBrewer)
create_spatial_feature <- function(pdata,feat.name,
                                   xloc.name = "sdimx",yloc.name = "sdimy",point.size = 0.5,
                                   color.scale = scale_colour_gradient(low = "grey",high = "red"))
{
  require(ggplot2)
  if (is.numeric(pdata[[feat.name]]))
  {
    if (is.null(color.scale))
    {
      color.scale = scale_colour_gradient(low = "grey",high = "red")
    }
    hobj = ggplot() + 
      geom_point(data = pdata,aes(x = .data[[xloc.name]],y = .data[[yloc.name]],colour = .data[[feat.name]]),size = point.size) + 
      color.scale + 
      theme_classic() + 
      theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())
    
  }
  if (is.character(pdata[[feat.name]]) | is.factor(pdata[[feat.name]]))
  {
    if (is.null(color.scale))
    {
      color.scale = scale_colour_brewer(palette = "Dark2")
    }
    hobj = ggplot() + 
      geom_point(data = pdata,aes(x = .data[[xloc.name]],y = .data[[yloc.name]],colour = .data[[feat.name]]),size = point.size) + 
      color.scale + 
      theme_classic() + 
      theme(axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())
  }
  return(hobj)
}

##### set paths
root.dir = "/media/won-min/My Passport1/SingleCell/Cancer/Melanoma/Data_for_ZNF180_Manuscript";
setwd(root.dir)

# params for neighborhood analysis
ex.thresh = 1 # expression threshold to define voxels with valid gene expressions
knn = 50 # number of k-nearest neighbors

##### Get reference scRNA-seq for InstaPrism run
# built-in reference from InstaPrism
refPhi_obj = InstaPrism_reference("SKCM")

## load pre-processed visium data
filenames <- list.files(path = root.dir,pattern = "processed_giotto.RDS",recursive = T,full.names = T)
names(filenames) = sapply(strsplit(filenames,"/"),function(x) x[10])

# assign files to save cell type deconvolution outputs
deconvfiles = gsub("processed_giotto.RDS","processed_giotto.InstaPrism_wt_SKCM.RDS",filenames)

## Run cell type identification + enrichments/depletion analysis in ZNF180 tumor neighborhood
stat.data = fet.data = data.frame()
ii = file.exists(deconvfiles)
filenames = filenames[ii]
deconvfiles = deconvfiles[ii]
print(deconvfiles)
spat.plots = vector("list",length(filenames))
names(spat.plots) = names(filenames)
for (i in 1:length(filenames))
{
  ## load spatial transcriptome
  cat(paste0("processing:",filenames[i],"\n"))
  visium = readRDS(filenames[i])
  visium = createSpatialNetwork(visium,method = "kNN",k = knn) # create spatial network using the kNN parameter. This will later define the cellular neighborhood as directly connected voxels. 
  
  # grab data
  if (TRUE)
  {
    cnt = visium@expression$cell$rna$raw@exprMat
    norm.dat = visium@expression$cell$rna$normalized@exprMat
    cell.feats = pDataDT(visium)
    spatial_locations = as.data.frame(visium@spatial_locs$cell$raw@coordinates)
  }
  
  # use InstaPrism to get cell type abundance 
  deconvfile = deconvfiles[i]
  if (!file.exists(deconvfile))
  {
    # Run InstaPrism
    deconv_res = InstaPrism(bulk_Expr = as.matrix(cnt), refPhi_cs = refPhi_obj,n.core = 10)
    saveRDS(deconv_res,file = deconvfile)
  }else{
    deconv_res = readRDS(file = deconvfile)
  }
  estimated_frac = t(deconv_res@Post.ini.cs@theta)
  estimated_frac = estimated_frac[match(colnames(norm.dat),rownames(estimated_frac)),];
  
  # add Melanoma cells
  mel.i = grep("Melanoma cells",colnames(estimated_frac))
  estimated_frac = cbind(estimated_frac[,setdiff(1:ncol(estimated_frac),mel.i)],rowSums(estimated_frac[,mel.i]))
  colnames(estimated_frac)[ncol(estimated_frac)] = "Melanoma cells"
  colnames(estimated_frac) = c("Cancer associated fibroblasts" = "CAFs","B cells" = "B cells","CD8+ T cells" = "CD8+ T cells","CD4+ T cells" = "CD4+ T cells",
                               "regulatory T cells" = "Tregs","Natural killer cells" = "NK cells","Macrophages/Monocytes" = "Macrophages/Monocytes",
                               "Dendritic cells" = "DCs", "Endothelial cells" = "Endothelial", "Melanoma cells" = "Melanoma")[colnames(estimated_frac)]
  
  celltype.binary = estimated_frac > 0.33 # binarized cell type presence with 1/3 threshold for Slide-Seq V2 platforms
  
  # identify cell type call: if multiple cell types are present in a voxel, concatenate their names
  ct = apply(celltype.binary,1,function(x) {out = NA;if (any(x)) out = paste(names(x)[x],collapse = "/");return(out)})
  
  # Define ZNF180+ tumor locations
  if (TRUE & any(rownames(norm.dat) == "ZNF180"))
  {
    #ii = colnames(norm.dat)[(norm.dat["ZNF180",] > ex.thresh) & estimated_frac[,colnames(estimated_frac) %in% c("Tumor","Melanoma cells")] > 0.5]
    ii = colnames(norm.dat)[(norm.dat["ZNF180",] > ex.thresh)]
    
    # define neighborhood
    nb.list = findNetworkNeighbors(gobject = visium,spatial_network_name = "kNN_network",source_cell_ids = ii)
    nb.list = nb.list[match(colnames(norm.dat),nb.list$cell_ID),]
    pdat = cbind.data.frame(spatial_locations,nb.list[match(nb.list[[1]],spatial_locations$cell_ID),-1]);
    pdat$celltype = ct[match(pdat$cell_ID,names(ct))]
    
    # refine neighbor call
    vec = rep(NA,nrow(pdat))
    vec[pdat$nb_cells %in% c("both","source")] = "ZNF180+ Tumor"
    vec[pdat$nb_cells == "neighbor"] = "Neighbor"
    #vec[estimated_frac[,"CD4+ T cells"] > 0.33 ] = "CD4+ T cells"
    
    if (any(rownames(norm.dat) == "TIGIT")) vec[estimated_frac[,"CD4+ T cells"] > 0.33 & norm.dat["TIGIT",] > ex.thresh] = "TIGIT+ CD4+ T cells"
    if (any(rownames(norm.dat) == "CD226")) vec[estimated_frac[,"CD4+ T cells"] > 0.33 & norm.dat["CD226",] > ex.thresh] = "CD226+ CD4+ T cells"
    
    vec = factor(vec,levels = c("ZNF180+ Tumor","Neighbor","CD4+ T cells","TIGIT+ CD4+ T cells","CD226+ CD4+ T cells"))
    pdat$neighbor.call = vec
    cols = c("ZNF180+ Tumor" = "red","Neighbor" = "goldenrod2","others" = "gray90","CD4+ T cells" = "black","TIGIT+ CD4+ T cells"= "magenta","CD226+ CD4+ T cells"= "blue")
    
    # make neighbor plot
    pdata = pdat;
    feat.name = "neighbor.call";
    xloc.name = "sdimx";
    yloc.name = "sdimy";
    point.size = 0.3;
    
    p.gene = ggplot() + 
      geom_point(data = pdata,aes(x = .data[[xloc.name]],y = .data[[yloc.name]],colour = .data[[feat.name]],alpha = .data[[feat.name]]),size = point.size) + 
      labs(title = names(filenames)[i]) + 
      scale_colour_manual(values = cols,na.value = "gray90",limits = c("ZNF180+ Tumor","Neighbor","CD4+ T cells","TIGIT+ CD4+ T cells","CD226+ CD4+ T cells")) + 
      scale_alpha_manual(values = c("ZNF180+ Tumor" = 0.85,"Neighbor" = 0.75,"others" = 0.2,"CD4+ T cells" = 0.85,"TIGIT+ CD4+ T cells" = 0.85,"CD226+ CD4+ T cells" = 0.85),na.value = 0.05) + 
      theme_classic() + 
      guides(alpha = "none",colour = guide_legend(override.aes = list(size = 5))) + 
      theme(plot.title = element_text(hjust = 0.5,size = 15),axis.title = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.line = element_blank(),legend.title = element_blank(),
            legend.text = element_text(size = 14),legend.position = "bottom")
    
    spat.plots[[i]] = p.gene
  }

  ### check enrichment/depletion of cell types in neighborhoods
  if (TRUE & any(rownames(norm.dat) == "ZNF180"))
  {
    source("scripts/R_functions/enrichment_functions.v3.R")
    
    # ZNF180 neighborhood set
    nb.set = subset(nb.list,nb_cells != "others" & !is.na(nb_cells))$cell_ID
    ct.set = apply(estimated_frac,2,function(x) names(x)[x > 0.3])
    sapply(ct.set,length)
    ct.set = ct.set[sapply(ct.set,length) > 5]
    fet.res = perform.AllPairs.FET(geneSets1 = list("ZNF180 Tumor NB" = nb.set),geneSets2 = ct.set,background = colnames(norm.dat),or = 1,
                                   alternative = "two.sided",adjust.FET.pvalue = T,do.multicore = F,n.cores = NULL)
    fet.res$data.id = names(filenames)[i]
    fet.res$corrected.FET.pvalue = p.adjust(fet.res$FET_pvalue,"BH")
    fet.data = rbind.data.frame(fet.data,fet.res)
  }
  
 
}

#### save output plots and tables
if (TRUE)
{
  fet.data$kNN = knn
  saveRDS(fet.data,file = paste0("ZNF180_neighborhood.enrichment_knn",knn,".RDS"))
  saveRDS(spat.plots,file = paste0("ZNF180_neighborhood.enrichment_knn",knn,"_ggplots.RDS"))
}else{
  fet.data = readRDS(file = paste0("ZNF180_neighborhood.enrichment_knn",knn,".RDS"))
  spat.plots = readRDS(file = paste0("ZNF180_neighborhood.enrichment_knn",knn,"_ggplots.RDS"))
}

#### Make examplary plots for manuscript
library(cowplot)
library(ggpubr)
plst = spat.plots[c("ECM08","MBM08","MBM06","MBM11_rep2")]
plst = lapply(plst,function(x) x + theme(legend.text = element_text(size = 18),plot.title = element_text(size = 20)))
plst[[1]] = plst[[1]] + theme(legend.text = element_text(size = 18),plot.title = element_text(size = 20))

lgd = as_ggplot(get_legend(plst[[1]]))
plst = lapply(plst,function(x) x + guides(colour = "none"))

png(file = "/media/won-min/My Passport1/SingleCell/Cancer/Melanoma/sc_Biermann_et_al_2022/spatial/ZNF180_neighborhood.png",res = 500,height= 4500,width = 3700)
plot_grid(plotlist = plst,ncol = 2)
#plot_grid(plot_grid(plotlist = plst,ncol = 2),lgd,ncol = 1,rel_heights = c(0.9,0.1))
dev.off()

##### 
#### Create statistics plots
library(ggplot2)
library(ggbeeswarm)
library(reshape2)
pdata = fet.data
pdata$odds.ratio[pdata$odds.ratio > 2] = 2
pdata$tumor.type = gsub("[0-9](.*)$","",pdata$data.id)

# create boxplot
pdata = fet.data
pdata$tumor.type = gsub("[0-9](.*)$","",pdata$data.id)
pdata$set2_Name = gsub("/","/\n",pdata$set2_Name)
boxobj = ggplot() + 
  #geom_quasirandom(data = pdata,aes(x = set2_Name,y = odds.ratio,colour =tumor.type,group = tumor.type),alpha = 0.5,position = position_jitterdodge()) + 
  #geom_boxplot(data = pdata,aes(x = set2_Name,y = odds.ratio,fill = tumor.type,colour = tumor.type),alpha = 0.5,outlier.colour = NA) + 
  geom_quasirandom(data = pdata,aes(x = set2_Name,y = odds.ratio),alpha = 0.5,position = position_jitterdodge()) + 
  geom_boxplot(data = pdata,aes(x = set2_Name,y = odds.ratio),alpha = 0.5,outlier.colour = NA) + 
  #scale_x_discrete(limits = ctick) +
  scale_y_continuous(limits = c(0,3)) + 
  geom_hline(yintercept = 1,colour = "red",linetype = "dashed") + 
  theme_classic() + 
  #guides(colour = guide_legend(title = "Metastasis\nSite"),fill = "none") + 
  theme(axis.text.x = element_text(size = 15,angle = 45,vjust = 1,hjust = 1),axis.title.y = element_blank(),axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 16),
        legend.title = element_text(size = 16),legend.text = element_text(size = 13),legend.position = "bottom",legend.direction = "vertical") + 
  coord_flip() + 
  labs(y = "Odds Ratio") 

png(file = "ZNF180Tumor_NB_vs_CellTypes.OR.png",res = 450,width = 3000,height = 2000)
print(boxobj)
dev.off()

pdf(file = "Spatial_Summary.pdf",width = 9,height = 6)
plot_grid(boxobj,plot_grid(plotlist = plst,ncol = 2),ncol = 2,rel_widths = c(0.4,0.6))
dev.off()

# spatial plots
library(cowplot)
png(file = "CD4_Tcell_vs_ZNF180_Tumor.png",res = 450,width = 5000,height = 4000)
plot_grid(plotlist = spat.plots[subset(fet.data,set2_Name == "CD4+ T cells")$data.id],ncol = 2)
dev.off()