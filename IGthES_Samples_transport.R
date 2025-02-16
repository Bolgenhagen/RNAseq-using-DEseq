library(rhdf5)
library(tximportData)
library(tximport)
library(readr)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(PCAtools)
library(EnhancedVolcano)
library("pheatmap")
library(openxlsx)

#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("PCAtools")

#is_loaded <- library("tximport", logical.return = TRUE)
#if (is_loaded) {
#  print("Library loaded successfully!")
#} else {
#  print("Library failed to load.")
#}


dir = ""
setwd (dir)

samples  <- list.files()
numeric_part <- as.numeric(gsub("\\D", "", samples)) 
ordered_samples <- samples[order(numeric_part)]

files <- file.path(dir, ordered_samples, "abundance.h5")
names(files) <- ordered_samples

ensembl <- useEnsembl(biomart = "genes")
searchDatasets(mart = ensembl, pattern = "dlabrax2021")

my_ensembl <- useEnsembl(biomart = "genes", dataset = "dlabrax_gene_ensembl")

attributes = listAttributes(my_ensembl)
head(attributes)

transcript2gene <- getBM(attributes=c("ensembl_transcript_id_version", "ensembl_gene_id"), mart = my_ensembl, values="")


txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = transcript2gene)

phenotypes = data.frame(matrix(vector(), 16, 2, dimnames=list(c(), c("Time", "Treatment"))),stringsAsFactors=F)
row.names(phenotypes) <- names(txi.kallisto.tsv$infReps)

phenotypes$Time <- c(rep(0,4), rep(4,4), rep(24, 8))
phenotypes$Treatment <- c(rep("Without Stress",4), rep("With Stress",4), rep("Without Stress",4), rep("With Stress",4))

data <- DESeqDataSetFromTximport(txi = txi.kallisto.tsv, colData = phenotypes, design = ~ as.factor(Time) + as.factor(Treatment))


vsd <- varianceStabilizingTransformation(data)

plotPCA(vsd, intgroup = c("Time","Treatment"))

p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)

screeplot(p)   #Percentage explained by each PCA

biplot(p, legendPosition = 'top', lab = NULL, legendLabSize = 12, legendIconSize = 6, title = "Test", pointSize=5)

pairsplot(p, components = c(1:5), triangle = TRUE, trianglelabSize = 12, hline = 0, vline = 0, pointSize = 1, gridlines.major = FALSE, gridlines.minor = FALSE, title = 'Pairs plot', plotaxes = FALSE, margingaps = unit(c(-0.01, -0.01, -0.01, -0.01), 'cm'))

data$combined <- factor(paste0(data$Time, data$Treatment))

design(data) <- ~ combined

dds <- DESeq(data)

res1 <- results(dds, contrast=c("combined","4With Stress","0Without Stress"))
res2 <- results(dds, contrast=c("combined","24With Stress","0Without Stress"))



annotation <- getBM(attributes=c("ensembl_gene_id", "external_gene_name","description"), mart = my_ensembl, values="")
row.names(annotation) <- annotation$ensembl_gene_id
annotation <- annotation[-1]

res1_annotated <- merge(DataFrame(res1), annotation, all.x = TRUE, by = "row.names")
res1_annotated_filtered <- subset(res1_annotated[order(res1_annotated$padj),], padj < 0.05)

res2_annotated <- merge(DataFrame(res2), annotation, all.x = TRUE, by = "row.names")
res2_annotated_filtered <- subset(res2_annotated[order(res2_annotated$padj),], padj < 0.05)




write.table(as.data.frame(res1_annotated_filtered),file="DE-results1.txt", sep="\t", row.names = FALSE)
write.table(as.data.frame(res2_annotated_filtered),file="DE-results2.txt", sep="\t", row.names = FALSE)
)


EnhancedVolcano(res1_annotated_filtered, x="log2FoldChange", y="padj", lab = "")
EnhancedVolcano(res2_annotated, x="log2FoldChange", y="padj", lab = "")


top501 <- res1_annotated_filtered[1:50,1]

top501_expression <- assay(vsd[top501])

top501_expression  <- top501_expression - rowMeans(top501_expression)

pheatmap(top501_expression)


top502 <- res2_annotated_filtered[1:50,1]

top502_expression <- assay(vsd[top502])

top502_expression  <- top502_expression - rowMeans(top502_expression)

pheatmap(top502_expression)




annot <- as.data.frame(colData(vsd)[, c("Treatment","Time")])


gene_names1 <- res1_annotated_filtered[1:50, "Row.names"]
pheatmap(top501_expression, annotation_col = annot)

gene_names2 <- res2_annotated_filtered[1:50, "Row.names"]
pheatmap(top502_expression, annotation_col = annot)


dds$Temperature <- as.factor(dds$Treatment)


dds$Time <- as.factor(dds$Time)

####### Aqui e para se eu quiser somente um gene em especifico mostra a expressao dele ao longo do experimento 
plotCounts(dds, gene= res1_annotated_filtered[1,1] , intgroup=c("Treatment","Time"))  

d <- plotCounts(dds, gene=res1_annotated_filtered[1,1], intgroup=c("Treatment","Time"), returnData = TRUE)
ggplot(d, aes(x=Treatment, y=count, col=Time)) + geom_point(position=position_jitter(w=0.1,h=0), size=1.5) + scale_y_log10() + theme_bw()

#### BASED ON  https://github.com/Roslin-Aquaculture/RNA-Seq-kallisto?tab=readme-ov-file#4-differential-expression-using-DESeq2
