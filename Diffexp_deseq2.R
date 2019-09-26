
#----------------------------ReadME use DESeq2----------------------------------------
## try http:// if https:// URLs are not supported

# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")


# cts -> is the matrix counts
# coldata -> is the sample information

# NOTE: "cts" in the first column must have "gene_id" and not be empty!
# NOTE2: When import "coldata", you must specified separation by \t

# pasCts <- read.csv("a_n0_m1_first_run.csv",sep="\t",row.names="gene_id")

#---------------------------Run DESeq2 ------------------------------------------------

library("DESeq2")
args = commandArgs(trailingOnly=TRUE)
cts <- as.matrix(read.csv(args[1],row.names="gene_id", header = T,sep = "\t"))
coldata <- read.csv(args[2], header = T, row.names = 1)

# With the count matrix, cts, and the sample information, coldata, 
# we can construct a DESeqDataSet:

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~lifestage)



# CONDITIONS FROM SAMPLE CONDITION FILE
dds$lifestage <- factor(dds$lifestage, levels = c("amastigote", "metacyclic_promatigote", "procyclic_promastigote"))


#run Deseq
dds <- DESeq(dds)
#---------------------------Results all vs all -----------------------------------------------
#obtain Table
res <- results(dds)
write.table(as.data.frame(res), 
            file="LbraM2903_DESeq2.csv",
            sep = "\t")
head(results(dds, tidy=TRUE))


#Summary of differential gene expression
summary(res) 

#Sort summary list by p-value

res <- res[order(res$padj,res$log2FoldChange),]
q=head(res, n=500
     )
write.table(as.data.frame(q), 
            file="LbraM2903_2_DESeq2.csv",
            sep = "\t")
head(results(dds, tidy=TRUE))



#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log2(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log2(pvalue), pch=20, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log2(pvalue), pch=20, col="red"))

#PCA plot all conditions 
library(ggplot2)
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="lifestage")
vsdata
# --------------------------COMPARATIVE VS lifestage----------------------------
res2 <- results(dds, contrast=c("lifestage", "amastigote", "metacyclic_promatigote"))
write.table(as.data.frame(res2), 
            file="lifestage_comparison_LbraM2903_amastigote_metacyclic.csv",
            sep = "\t")

res3 <- results(dds, contrast=c("lifestage","amastigote","procyclic_promastigote"))
write.csv(as.data.frame(res3), 
          file="lifestage_comparison_LbraM2903_amastigote_procyclic.csv")

res4 <- results(dds, contrast=c("lifestage","metacyclic_promatigote","procyclic_promastigote"))
write.csv(as.data.frame(res4), 
          file="lifestage_comparison_LbraM2903_metacyclic_procyclic.csv.csv")


#---------------------------Normalized count -----------------------------------------------


dds1 <- estimateSizeFactors(res2)
dds2<-counts(dds1, normalized=TRUE)
write.csv(as.data.frame(dds2), file = "norm_reads_LbraM2903.csv")

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)



#-----------------------------PLOT Heatmap of DE genes-------------------------------------------------


#metadata
my_sample_col2<- data.frame(rownames(coldata), coldata$lifestage)
colnames(my_sample_col2) <- c("run","lifestyle")
rownames(my_sample_col2)<- my_sample_col2$run
my_sample_col2$run<- NULL
my_colour = list(lifestyle=c(amastigote="#1a5276",metacyclic_promatigote="#ff7f00", 
                             procyclic_promastigote ="#d4ac0d"))

#heatmap
library(pheatmap)

resSort <- res[order(res$padj, res$log2FoldChange),]
head(resSort)


topgenes <- head(rownames(subset(res, padj < 0.05 | log2FoldChange < -1.5 |
                                   log2FoldChange > 1.5)),500)

write.table(as.data.frame(topgenes), 
            file="LbraM2903_topgenes.csv",
            sep = "\t")

rld <- rlog(dds, blind=FALSE)
mat <- assay(rld)[topgenes,]
mat <- mat - rowMeans(mat)
pheatmap(mat, 
         annotation_colors = my_colour,
         annotation_col=my_sample_col2,
         fontsize = 8)




par(mfrow=c(2,3))

plotCounts(dds, gene="LbraM2903_20_1855511_1855532", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_20_1563957_1563978", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_20_1890086_1890107", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_18_87666_87686", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_23_304421_304441", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_26_687143_687164", intgroup="lifestage")


par(mfrow=c(2,3))

plotCounts(dds, gene="LbraM2903_36_2664176_2664257", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_34_94566_94639", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_31_1249411_1249434", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_14_343046_343079", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_9_394024_394045", intgroup="lifestage")
plotCounts(dds, gene="LbraM2903_28_963125_963148", intgroup="lifestage")



# volcano plot 

res <- read.table("LbraM2903_DESeq2.csv", header=TRUE, sep = "\t", na.strings = )
head(res)

# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>2), textxy(log2FoldChange, -log10(pvalue), labs=gene_ID, cex=.8))
