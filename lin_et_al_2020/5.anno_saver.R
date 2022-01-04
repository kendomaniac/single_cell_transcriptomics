WD = "/data/single_cell_transcriptomics/test"
INPUT = "matrix.csv"
SEPARATOR = ","
GTF = "Homo_sapiens.GRCh38.105.gtf"
SCRATCH = "/scratch"
setwd(WD)
library(rCASC)

scannobyGtf(
  group = "docker", 
  file=paste(getwd(), paste("saver", INPUT, 
             sep="_"), sep="/"), 
  gtf.name = GTF, 
  biotype = "protein_coding", 
  mt = FALSE, 
  ribo.proteins = FALSE, 
  umiXgene = 3, 
  riboStart.percentage = 0, 
  riboEnd.percentage = 100, 
  mitoStart.percentage = 0, 
  mitoEnd.percentage = 100, 
  thresholdGenes = 1)

raw <- read.table(paste("filtered_annotated", INPUT, 
                  sep="_"), 
                  sep=SEPARATOR, header=T, row.names=1)
saver <- read.table(paste("annotated_saver", INPUT, 
                    sep="_"), 
                  sep=SEPARATOR, header=T, row.names=1)
saver_ribomito <- saver[,which(names(saver)%in%names(raw))]
write.table(saver_ribomito, paste("filtered_annotated_saver_ribomito", INPUT, sep="_"), sep=",", col.names=NA)

