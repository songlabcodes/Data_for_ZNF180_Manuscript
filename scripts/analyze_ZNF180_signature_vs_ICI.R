#### create figure for 
rm(list = ls())

library(readxl)
library(cowplot)
library(Seurat)

root.dir <- getwd();setwd(root.dir)
encodefile = "ZNF180_peaks.anno.txt" # ENCODE ZNF180 TF binding data
znf180.sigfile = "Song_et_al_2021/41467_2021_21457_MOESM9_ESM.xlsx" # signatures from publication
liu.folder = "Liu_et_al_2019" # folder to hold Liu et al. 2019 data
riaz.folder = "Riaz_et_al_2016" # folder to hold Riaz et al. 2016 data
########## Summarize ENCODE data
# make pie chart in HepG2 cells
if (TRUE)
{
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(reshape)
  
  encode.peaks <- read.delim(file = encodefile,sep = "\t",header = TRUE,stringsAsFactors = FALSE)
  encode.peaks$broad.annotation = gsub("Exon \\((.*)$","Exon",gsub("Intron \\((.*)$","Intron",encode.peaks$annotation))
  
  pdat = melt(table(encode.peaks$broad.annotation));colnames(pdat) = c("Type","Count")
  pdat$Count = pdat$Count/sum(pdat$Count)
  
  pdat <- pdat %>% 
    arrange(desc(Type)) %>%
    mutate(prop = Count / sum(pdat$Count)) %>%
    mutate(ypos = cumsum(prop)- 0.5*prop )
  
  pie.obj = ggplot(data = pdat, aes(x="", y=Count, fill=Type)) +
    geom_bar(stat="identity", width=1,colour = "white") +
    scale_y_continuous(labels = scales::percent) + 
    coord_polar("y", start=0) + theme_classic() + 
    theme(axis.title = element_blank(),axis.line = element_blank(),axis.text = element_text(size = 15),
          legend.text = element_text(size = 14),legend.title = element_blank(),legend.position = "bottom")+
    scale_fill_brewer(palette="Set1") + 
    guides(fill = guide_legend(ncol = 2))
  
  print(pie.obj)
}


# siZNF180 signature: SKMEL147. Directly import the signature from Song et al. 2021
if (TRUE)
{
  znf180.sig.table = read_excel(path = znf180.sigfile,sheet = 1,skip = 2)
  df = subset(znf180.sig.table,adj.P.Val < 0.05 & abs(logFC) > log2(1.5) & target.gene == "ZNF180")
  siRNA.sigs = split(df$ID,df$signature.id)
  names(siRNA.sigs) = paste0("SKMEL147.si",names(siRNA.sigs))
}

######## Check if ZNF180 signatures are predicted of response 
if (TRUE)
{
  source("scripts/R_functions/survival_multivariate_analysis.R")
  
  # load Liu et al. 2019 data
  if (TRUE)
  {
    # load processed data by the published work
    data.df = read.delim(file= paste0(liu.folder,"/41591_2019_654_MOESM3_ESM.txt"),sep = "\t",header = T,stringsAsFactors = F)
    data.mat= as.matrix(data.df[,-1])
    rownames(data.mat) = data.df[[1]]
    data.mat = t(data.mat)
    data.mat = log2(data.mat+1)
    rownames(data.mat) = gsub("HLA\\.","HLA-",rownames(data.mat))
    
    library(readxl)
    # load clinical follow up data
    cif.df = read_excel(path = paste0(liu.folder,"/41591_2019_654_MOESM4_ESM.xlsx"),sheet = 1,skip = 2)
    cif.df = cif.df[match(colnames(data.mat),cif.df[[1]]),]
  }
  
  # get intersection of ENCODE TF binding and Song et al. 2021 signatures
  znf180.sig = list(ZNF180.Regulated = intersect(siRNA.sigs$`SKMEL147.siZNF180-NTC_DN`,subset(encode.peaks,broad.annotation == "Promoter (<=1kb)")$SYMBOL))
  
  # run GSVA
  library(GSVA)
  gsva.score = gsva(expr = data.mat,gset.idx.list = znf180.sig,method = "ssgsea")
  cif.wt.gsva = cif.df
  cif.wt.gsva$ZNF180.Regulated = gsva.score[,cif.df[[1]]]
  
  km.out = run_kmplot(sdf = cif.wt.gsva,gene.id = "ZNF180.Regulated",time.col = "OS",event.col = "dead",qtiles = 0.5)
  km.obj = km.out$kmplot + labs(y = "Survival Probability") + 
    guides(colour = guide_legend("ZNF180 Regulation\nby ssGSEA")) + 
    theme_classic() + theme(axis.text = element_text(size = 14),axis.title = element_text(size = 17),
                            plot.title = element_text(size = 19,hjust = 0.5),legend.position = "bottom",legend.title.position = "left",
                            legend.title = element_text(size = 15),legend.text = element_text(size = 13))
  print(km.obj)
  #### load Riaz et al. 2016 data
  if (TRUE)
  {
    library(GSVA)
    library(grid)
    dge.obj = readRDS(file = paste0(riaz.folder,"/DGEList.RDS")) # load pre-processed data by edgeR workflow to CPM + TMM normalization
    coldat = dge.obj$samples
    coldat$treatment = coldat$visit..pre.or.on.treatment..ch1
    coldat = subset(coldat,response.ch1 != "UNK")
    #coldat = subset(coldat,treatment == "Pre")
    coldat$response.ch1 = factor(coldat$response.ch1)
    coldat$exp.design.group = paste0(coldat$treatment,"_",coldat$response.ch1)
    
    dge.obj = calcNormFactors(dge.obj)
    data.mat = cpm(dge.obj,log = T)
    
    gsva.score = gsva(expr = data.mat,znf180.sig,method = "ssgsea")
    
    ### Perform analysis 
    #df =subset(coldat,treatment == "On")
    df = coldat
    m = rbind(gsva.score[,match(rownames(df),colnames(gsva.score))])
    rownames(m) = names(znf180.sig)
    
    case.ii = which(df$response.ch1 == "PD")
    ctrl.ii = which(df$response.ch1 %in% c("PRCR"))
    stat.res = apply(rbind(m),1,function(x) {out = wilcox.test(x = x[case.ii],y = x[ctrl.ii]);c("p.value" = out$p.value,mean.case = mean(x[case.ii]),mean.ctrl = mean(x[ctrl.ii]))})
    
    # show boxplots
    library(ggbeeswarm)
    pdata = cbind.data.frame(df,as.data.frame(t(m)))
    pdata = subset(pdata,response.ch1 %in% c("PD","PRCR"))
    
    stat.chr = paste0("Wiconxon P=",format(signif(stat.res[1,1],3),scientific = T))
    pobj = ggplot(data = pdata) + geom_quasirandom(aes(x = response.ch1,y = ZNF180.Regulated,colour = treatment)) + 
      geom_boxplot(aes(x = response.ch1,y = ZNF180.Regulated),alpha = 0.2,width = 0.3) + 
      theme_classic() + 
      guides(colour = "none") + 
      labs(y = "GSVA Score") + 
      theme(axis.text.x = element_text(size = 18),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 18),
            axis.title.y = element_text(size = 17),
            legend.position = "right",legend.direction = "vertical",
            legend.title = element_text(size = 14),legend.text = element_text(size = 13))
    
    # Create a text
    grob <- grobTree(textGrob(stat.chr, x=0.5,  y=0.965, hjust=0.5,
                              gp=gpar(col="black", fontsize=19)))
    
    pobj = pobj + annotation_custom(grob)
    
  }
  
  prog.obj = plot_grid(km.obj,pobj,ncol = 1,rel_heights = c(0.6,0.4))
}

print(prog.obj)