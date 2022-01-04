WD = "/data/single_cell_transcriptomics/test"
INPUT = "matrix.csv"
SEPARATOR = ","
GTF = "Homo_sapiens.GRCh38.105.gtf"
SCRATCH = "/scratch"
setwd(WD)
library(rCASC)
scannobyGtf(
  group = "docker",
  file = paste(getwd(), INPUT, sep="/"),
  gtf.name = GTF,
  biotype = "protein_coding",
  mt = FALSE,
  ribo.proteins = FALSE,
  umiXgene = 3,
  riboStart.percentage = 1,
  riboEnd.percentage = 50,
  mitoStart.percentage = 1,
  mitoEnd.percentage = 10,
  thresholdGenes = 101
)

