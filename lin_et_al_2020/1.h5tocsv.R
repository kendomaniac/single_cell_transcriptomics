library(rCASC)
WD = "/data/single_cell_transcriptomics/test"
INPUT = "matrix.mtx.gz"
SEPARATOR = ","
setwd(WD)
h5tocsv(
  group = "docker",
  file = paste(getwd(), INPUT, sep="/"),
  type = "10xgenomics",
  version = "5"
)

