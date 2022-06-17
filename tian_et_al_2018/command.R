data <- read.table("annotated_GSM3618014_gene_count.csv", sep=",", header=T, row.names=1)
clusters <- read.table("RNA-5c_clustering.output.csv", sep=",", header=T, row.names=1)
names(clusters)
identical(rownames(clusters), names(data))
common <- intersect(names(data), rownames(clusters))
data <- data[, which(names(data)%in%common)]
identical(rownames(clusters), names(data))
names(data) <- clusters$annotated
write.table(data, "RNA-5c.csv", sep=",", col.names=NA)

### saver and saver CPM are in 10.6084/m9.figshare.20087498 

library(SAVER)
WD = getwd()
INPUT = "RNA-5c.csv"
SEPARATOR = ","
THREADS = 12 # the number of cores dedicated to run such analysis
setwd(WD)
raw.data <- read.table(INPUT, sep=SEPARATOR, header = T, row.names=1)
dataset <- as.matrix(raw.data)
dataset.saver <- saver(dataset, ncores = THREADS, estimates.only = TRUE)
write.table(dataset.saver, paste("saver", INPUT, sep="_"), sep=SEPARATOR, col.names=NA)

# preprocessing
ann <- read.table("annotated_GSM3618014_gene_count.csv", sep=",", header=T, row.names=1)
symbol <- sapply(strsplit(rownames(ann), ":"), function(x)x[2])
id <- sapply(strsplit(rownames(ann), ":"), function(x)x[1])
df <- data.frame(symbol, id, stringsAsFactors=F)
rownames(df) <- (rownames(ann))
df <- df[which(df$id %in% rownames(saver)),]
saver <- read.table("saver_RNA-5c.csv", sep=",", header=T, row.names=1)
saver <- saver[order(rownames(saver)),]
df <- df[order(df$id),]
identical(rownames(saver), df$id)
saver <- saver[!duplicated(df$symbol),]
df <- df[!duplicated(df$symbol),]
rownames(saver) <- df$symbol
write.table(saver, "symbol_annotated_saver_RNA-5c.csv", sep=",", col.names=NA)

counts2cpm <- function(file, sep = ","){
        tmp <- read.table(file, sep=sep, header=T, row.names=1)
        col.sum <- apply(tmp, 2, sum)
        tmp1 <- t(tmp)/col.sum
        tmp1 <- t(tmp1)
        tmp1 <- tmp1 * 1000000
        write.table(tmp1, "symbol_annotated_saver_RNA-5c_cpm.csv", sep=",", 
 col.names=NA)
}

counts2cpm(file="symbol_annotated_saver_RNA-5c.csv", sep=",")
             
#saver and saver cpm data are available in 10.6084/m9.figshare.20087498

# sca run version 1 or 2
library(rCASC)
autoencoder4clustering(
  group = "docker",
  scratch.folder="/scratch",
  file=paste(getwd(), "symbol_annotated_saver_RNA-5c_cpm.csv", sep="/"),
  separator=",",
  bias="cytoBands",
  permutation=40,
  nEpochs=200,
  patiencePercentage = 5,
  seed = 1111,
  projectName="CYTOcpu",
  bN = "NULL",
  lr = 0.01,
  beta_1 = 0.9,
  beta_2 = 0.999,
  epsilon = 1e-08,
  decay = 0,
  loss = "mean_squared_error",
  regularization = 10,
  version=1
)
             
# post processing (~/Results/CYTO/permutation)
files=list.files(getwd())
files-files[grep("denseSpace", files)]
temp=read.table(files[1],sep=",",header=TRUE,row.names=1)
for(i in files[-1]){
  temp=temp+read.table(i,header=TRUE,row.names=1,sep=",")
}
write.table(temp,"total.csv",col.names=NA,sep=",")
file=paste(getwd(),"total.csv",sep="/")
library(rCASC)
seuratBootstrap(group="docker", 
                scratch.folder="/scratch", 
                file=file, 
                nPerm=40, 
                permAtTime=5, 
                percent=10, 
                separator=",", 
                logTen=0, 
                pcaDimensions=20, 
                seed=1111,
                sparse=FALSE,
                format="NULL",
                resolution=0.2
)
             




