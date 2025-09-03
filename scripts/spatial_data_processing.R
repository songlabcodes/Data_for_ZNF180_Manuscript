rm(list = ls())

library(Giotto)
library(GiottoData)
require(Matrix)
library(readxl)

###########
root.dir <- "/media/won-min/My Passport1/SingleCell/Cancer/Melanoma/Data_for_ZNF180_Manuscript";
setwd(root.dir)
dir.create("Results")
out.dir <- "Results/Biermann_et_al_2022";dir.create(out.dir)

### Get list of input data 
file.data = readRDS(file = "Biermann_et_al_2022/file_annotation_wt_clinical.RDS")
file.data$coord.file = gsub("^(.*)/","Biermann_et_al_2022/SlideSeqV2/raw_data/",file.data$coord.file)
file.data$count.file = gsub("^(.*)/","Biermann_et_al_2022/SlideSeqV2/raw_data/",file.data$count.file)
########################## Data processing per sample
for (i in 1:nrow(file.data))
{
  
  cat(paste0("processing:",file.data$sample.name[i],"\n"))
  results_directory = paste0(out.dir,"/",file.data$sample.name[i])
  
  ### Set Giotto instructions
  my_python_path = NULL 
  instrs = createGiottoInstructions(save_plot = TRUE,
                                    show_plot = FALSE,
                                    save_dir = results_directory,
                                    python_path = my_python_path)
  
  ### load coordinates
  coord.dat = read.delim(file = gzfile(file.data$coord.file[i]),sep = ",",header = TRUE,stringsAsFactors = FALSE)[,-1]
  rw.nm = coord.dat[[1]]
  coord.dat = coord.dat[,-1];rownames(coord.dat) = rw.nm;rm(rw.nm)
  
  ### load count matrix
  expr.dat = read.delim(file = gzfile(file.data$count.file[i]),sep = ",",header = TRUE,stringsAsFactors = FALSE);
  rw.nm = as.character(expr.dat[[1]])
  expr.dat = expr.dat[,-1]
  rownames(expr.dat) = rw.nm
  rm(rw.nm)
  colnames(expr.dat) = gsub("\\.","-",colnames(expr.dat))
  expr.dat = Matrix(as.matrix(expr.dat))
  
  coord.dat = coord.dat[match(colnames(expr.dat),rownames(coord.dat)),]
  
  ### make meta data
  meta = data.frame(cell_ID = colnames(expr.dat),as.data.frame(file.data[rep(i,nrow(coord.dat)),c("patient id:ch1","mutation status:ch1","metastatic site:ch1")]))
  colnames(meta)[2] = "patient.ID"
  meta$sample.id = rep(file.data$sample.name[i],nrow(meta))
  
  ### create Giotto object\
  visium_dat <- createGiottoObject(expression = expr.dat, 
                                   spatial_locs = coord.dat, 
                                   instructions = instrs)
  
  # Add additional annotation if wanted
  visium_dat = addCellMetadata(visium_dat,
                               new_metadata = meta,
                               by_column = T,
                               column_cell_ID = 'cell_ID')
  
  ### filtering steps:
  library(cowplot)
  plot_grid(filterDistributions(visium_dat, detection = 'cells'),
            filterDistributions(visium_dat, detection = 'feats'))
  
  visium_dat <- filterGiotto(gobject = visium_dat,
                             expression_threshold = 1E-320,
                             gene_det_in_min_cells = 100,
                             min_det_genes_per_cell = 100)
  
  ### normalize
  # normalize to scale expression values of the Giotto object using the standard method, z-scoring feats over cells
  visium_dat <- normalizeGiotto(gobject = visium_dat,
                                norm_methods = 'standard',
                                scale_feats = TRUE,
                                scalefactor = 6000,
                                scale_order = 'first_feats', # Default, alternatively 'first_cells'
                                verbose = T)
  
  # Add gene & cell statistics to the Giotto object using the data normalized with the standard method
  visium_dat <- addStatistics(gobject = visium_dat, expression_values = 'normalized')
  head(pDataDT(visium_dat))
  head(fDataDT(visium_dat))
  
  
  mt_genes = grep("^MT-", x = fDataDT(visium_dat)$feat_ID, value = TRUE)
  ribo_genes = grep("^RPL|^RPS", x = fDataDT(visium_dat)$feat_ID, value = TRUE)
  
  visium_dat <- addFeatsPerc(visium_dat,
                             expression_values = 'normalized',
                             feats = mt_genes,
                             vector_name = "perc_mt")
  
  visium_dat <- addFeatsPerc(visium_dat,
                             expression_values = 'normalized',
                             feats = ribo_genes,
                             vector_name = "perc_ribo")
  
  # show MT and RIBO rates
  pobj = spatPlot2D(gobject = visium_dat, point_size = 1.5,cell_color = 'perc_mt', color_as_factor = F) + 
    spatPlot2D(gobject = visium_dat, point_size = 1.5,cell_color = 'perc_ribo', color_as_factor = F) + 
    spatPlot2D(gobject = visium_dat, point_size = 1.5,cell_color = 'nr_feats', color_as_factor = F) 
  
  png(file = paste0(results_directory,"/QC.metrics.png"),res =300,width = 8400,height = 2800)
  print(pobj)
  dev.off()
  # Since there are no known batch effects, the number of features detected per cell
  # will be regressed out so that covariates will not effect further analyses.
  visium_dat <- adjustGiottoMatrix(gobject = visium_dat,
                                   expression_values = c('normalized'),
                                   covariate_columns = c('nr_feats'))
  
  ##### Clustering
  # Calculate HVF using coefficient of variance within groups
  visium_dat <- calculateHVF(gobject = visium_dat, method = 'cov_groups')
  
  ## Select genes highly variable genes that fit specified statistics
  # These are both found within feature metadata
  feature_metadata = fDataDT(visium_dat)
  featgenes = feature_metadata[hvf == 'yes' & perc_cells > 1 & mean_expr_det > 0.5]$feat_ID
  
  ## run PCA on expression values (default)
  visium_dat <- runPCA(gobject = visium_dat, feats_to_use = featgenes, scale_unit = F, center = F)
  
  # add tSNE
  visium_dat <- runtSNE(visium_dat, dimensions_to_use = 1:15,check_duplicates = FALSE)
  
  # add UMAP
  visium_dat <- runUMAP(visium_dat, dimensions_to_use = 1:15)
  
  ## create a shared nearest neighbor network (sNN), where k is the number of k neighbors to use
  visium_dat <- createNearestNetwork(gobject = visium_dat, dimensions_to_use = 1:15, k = 15)
  
  ## Leiden clustering - increase the resolution to increase the number of clusters
  visium_dat <- doLeidenCluster(gobject = visium_dat,resolution = 0.1,
                                n_iterations = 1000,name = 'leiden_0.1_1000')
  
  #Plot UMAP post-clustering to visualize Leiden clusters
  umap.obj = plotUMAP(gobject = visium_dat,
                      cell_color = 'leiden_0.1_1000',
                      show_NN_network = F,
                      point_size = 1)
  
  ### Plot out clustering in 2d spatial dimension
  # Plot cell_color as a representation of the number of features/ cell ("nr_feats")
  spat.obj=spatPlot2D(gobject = visium_dat, point_size = 1,
                      cell_color = 'leiden_0.1_1000', color_as_factor = T)
  
  png(file = paste0(results_directory,"/Cluster_Results.png"),res =300,width = 5400,height = 2800)
  print(umap.obj + spat.obj)
  dev.off()
  
  saveRDS(visium_dat,file = paste0(results_directory,"/processed_giotto.RDS"))
  rm(visium_dat)
  gc()
}