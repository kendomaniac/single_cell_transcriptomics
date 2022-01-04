library(rCASC)
WD = "/data/single_cell_transcriptomics/test"
INPUT = "filtered_annotated_saver_ribomito_matrix.csv"
SEPARATOR = ","
SCRATCH = "/scratch"
setwd(WD)
seurat_ccycle(group = "docker", 
  scratch.folder = SCRATCH,
  file = paste(getwd(), INPUT, sep = "/"),     
  separator = SEPARATOR,
  seed=111)

