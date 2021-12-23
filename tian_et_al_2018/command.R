data <- read.table("annotated_GSM3618014_gene_count.csv", sep=",", header=T, row.names=1)
clusters <- read.table("RNA-5c_clustering.output.csv", sep=",", header=T, row.names=1)
names(clusters)
identical(rownames(clusters), names(data))
common <- intersect(names(data), rownames(clusters))
data <- data[, which(names(data)%in%common)]
identical(rownames(clusters), names(data))
names(data) <- clusters$annotated
write.table(data, "RNA-5c.csv", sep=",", col.names=NA)
