install.packages(c("BiocManager","pheatmap","RColorBrewer"))
BiocManager::install(c("airway","DESeq2"),force=TRUE)
library(airway)
library(DESeq2)
data(airway)

se <- airway
se

countdata <- assay(se)

coldata <- colData(se)

coldata

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design =
~ cell + dex)

counts(dds)

nrow(dds)

keep <- rowSums(counts(dds)) > 1

dds <- dds[keep, ]

nrow(dds) #29391

#####

#normalisation
rld <- rlog(dds, blind = FALSE)  
rld_countdata <- assay(rld)

plot(x = counts(dds)[,1], y = counts(dds)[,2], pch = 16)
plot(x = rld_countdata[,1], y = rld_countdata[,2], pch = 16)

##heatmaps

library(pheatmap)
library(RColorBrewer)

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(rld$cell, rld$dex, sep = "_")
colnames(sampleDistMatrix) <- paste(rld$cell, rld$dex, sep = "_")

colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pheatmap(mat = sampleDistMatrix, color = colors, clustering_distance_rows = 
         sampleDists, clustering_distance_cols = sampleDists)


##PCA

plotPCA(rld, intgroup = "cell", "dex")
