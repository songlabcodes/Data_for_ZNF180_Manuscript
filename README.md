### Overview
This github page contains codes to reproduce bioinformatics analysis in Song et al. 2025 to analyze melanoma -omics data. The input data are uploaded in Synapse (syn69530111). The data needs to be downloaded and extracted under the desired working folder. Once extracted, it is advised that "Results" folder is created under the current working directory to hold the outputs. Below is the list of codes: 

**analyze_ZNF180_signature_vs_ICI.R**: This code reproduces Figure 1B-E. 

**analyze_tumor_subclusters.R**: This reproduces Figure 1F-H. 

**PPI_analysis.R**: This reproduces PPI analysis in Figure 4A, C. 

**analyze_MMLines.R**: This reproduces Figure 4E to show ZNF180 co-expressions with AXL, FOSL1 in mesenchymal-like cells. 

**spatial_data_process.R**: This is the pre-processing code for SlideSeq-V2 data from Biermann et al. 2022 presented in Figure 5A, B. 

**analyze_ZNF180_tumor_spatial_transcriptome.R**: Using the pre-processed data from "spatial_data_process.R", this reproduces Figure 5A, B in the main text.

Make sure R function codes under /scripts/R_functions are available to call necessary functions. 