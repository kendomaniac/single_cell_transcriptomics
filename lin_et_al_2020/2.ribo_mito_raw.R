# Open an R session an type the following commands
WD = "/data/single_cell_transcriptomics/test"
INPUT = "matrix.csv"
SEPARATOR = ","
GTF = "Homo_sapiens.GRCh38.105.gtf"
SCRATCH = "/scratch"
setwd(WD)
library(rCASC)
mitoRiboUmi(
  group = "docker",
  scratch.folder = SCRATCH,
  file = paste(getwd(), INPUT, sep="/"),
  separator = ",",
  gtf.name = GTF,
  bio.type = "protein_coding",
  umiXgene = 3
)

