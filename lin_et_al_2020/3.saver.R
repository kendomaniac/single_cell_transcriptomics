library(SAVER)
WD = "/data/single_cell_transcriptomics/test"
INPUT = "matrix.csv"
SEPARATOR = ","
THREADS = 12 # the number of cores dedicated to run such analysis
setwd(WD)
raw.data <- read.table(INPUT, sep=SEPARATOR, header = T, row.names=1)
dataset <- as.matrix(raw.data)
dataset.saver <- saver(dataset, ncores = THREADS, estimates.only = TRUE)
write.table(dataset.saver, paste("saver", INPUT, sep="_"), sep=SEPARATOR, col.names=NA)

