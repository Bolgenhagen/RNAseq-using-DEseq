#### Load Required Libraries ####

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
library(grid)
library("showtext")
library(jsonlite)
library(gtools)
library(dplyr)
library(glmnet)
library(caret)
library(gtools)

#### Set Working Directory ####

dir = "C:/Users/felip/Documents/Mestrado/Tese/data_transport"
setwd(dir)

#### Extract Kallisto Mapping Statistics ####

# List all run_info.json files
files <- list.files(path = dir, pattern = "run_info.json", recursive = TRUE, full.names = TRUE)

# Extract number of processed reads
read_counts <- sapply(files, function(f) {
  info <- fromJSON(f)
  info$n_processed
})
# Print mean reads
mean_reads <- mean(read_counts)
cat("Average number of reads per sample:", round(mean_reads), "\n")
# Extract processed & pseudoaligned reads + calculate mapping %
mapping_stats <- do.call(rbind, lapply(files, function(f) {
  info <- fromJSON(f)
  data.frame(
    Sample = basename(dirname(f)),
    Processed_Reads = info$n_processed,
    Pseudoaligned_Reads = info$n_pseudoaligned,
    Mapping_Percent = round((info$n_pseudoaligned / info$n_processed) * 100, 2)
  )
}))
# Sort samples (natural order)
mapping_stats <- mapping_stats[mixedorder(mapping_stats$Sample), ]
rownames(mapping_stats) <- NULL
print(mapping_stats)
# Save mapping summary
write.csv(
  mapping_stats,
  "C:/Users/felip/Documents/Mestrado/Tese/Results/Transport/kallisto_mapping_summary.csv",
  row.names = FALSE
)

#### Prepare Files for tximport ####

samples  <- list.files()
numeric_part <- as.numeric(gsub("\\D", "", samples)) 
ordered_samples <- samples[order(numeric_part)]

files <- file.path(dir, ordered_samples, "abundance.h5")
names(files) <- ordered_samples

#### Prepare Biomart Annotation (Dlabrx dataset) ####

ensembl <- useEnsembl(biomart = "genes")
searchDatasets(mart = ensembl, pattern = "dlabrax2021")

my_ensembl <- useEnsembl(biomart = "genes", dataset = "dlabrax_gene_ensembl")

attributes = listAttributes(my_ensembl)
head(attributes)

transcript2gene <- getBM(
  attributes=c("ensembl_transcript_id_version", "ensembl_gene_id"),
  mart = my_ensembl, values=""
)

#### Import Quantifications with tximport ####

txi.kallisto.tsv <- tximport(files, type = "kallisto", tx2gene = transcript2gene)

#### Build Phenotype Table ####

phenotypes = data.frame(
  matrix(vector(), 16, 2,
         dimnames=list(c(), c("Time", "Treatment"))),
  stringsAsFactors=F)

row.names(phenotypes) <- names(txi.kallisto.tsv$infReps)

phenotypes$Time <- c(rep(0,4), rep(4,4), rep(24, 8))
phenotypes$Treatment <- c(rep("Without Stress",4),
                          rep("With Stress",4),
                          rep("Without Stress",4),
                          rep("With Stress",4))

#### DESeq2 Preprocessing ####

data <- DESeqDataSetFromTximport(
  txi = txi.kallisto.tsv,
  colData = phenotypes,
  design = ~ as.factor(Time) + as.factor(Treatment)
)

vsd <- varianceStabilizingTransformation(data)

#### PCA + Visualization (ggplot2 and PCAtools) ####

p1 <- plotPCA(vsd, intgroup = c("Time", "Treatment"))
p1 <- p1 + ggtitle("PCA Plot of Gene Expression") +
  labs(color = "Group", shape = "Group") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

ggsave("../Results/Transport/PCA_test.tiff", p1,
       width = 20, height = 16, dpi = 300, scale = 0.5)

p <- pca(assay(vsd), metadata = colData(vsd), removeVar = 0.1)

screeplot(p)
screeplot(p) + ggtitle("Scree Plot of PCA Analysis")

biplot(p, legendPosition = 'top', lab = NULL)

pairsplot(p, components = c(1:5))

#### Redefine Design for Combined Factor and Run DESeq ####

data$combined <- factor(paste0(data$Time, data$Treatment))
design(data) <- ~ combined

dds <- DESeq(data)

res1 <- results(dds, contrast=c("combined","4With Stress","0Without Stress"))
res2 <- results(dds, contrast=c("combined","24With Stress","24Without Stress"))

#### Annotate DE Results (Biomart) ####
  
annotation <- getBM(
  attributes=c("ensembl_gene_id","external_gene_name","description"),
  mart = my_ensembl, values=""
)

row.names(annotation) <- annotation$ensembl_gene_id
annotation <- annotation[-1]

res1_annotated <- merge(DataFrame(res1), annotation, all.x = TRUE, by = "row.names")
res1_annotated_filtered <- subset(res1_annotated[order(res1_annotated$padj),], padj < 0.05)

res2_annotated <- merge(DataFrame(res2), annotation, all.x = TRUE, by = "Row.names")
res2_annotated_filtered <- subset(res2_annotated[order(res2_annotated$padj),], padj < 0.05)

#### Volcano Plots ####

df1 = na.omit(res1_annotated)
df2 = na.omit(res2_annotated)

top_genes1 <- df1[df1$padj < 0.05 & abs(df1$log2FoldChange) > 1, ]
upregulated_genes1 <- sum(df1$log2FoldChange > 1 & df1$padj < 0.05)
downregulated_genes1 <- sum(df1$log2FoldChange < -1 & df1$padj < 0.05)

top_genes2 <- df2[df2$padj < 0.05 & abs(df2$log2FoldChange) > 1, ]
upregulated_genes2 <- sum(df2$log2FoldChange > 1 & df2$padj < 0.05)
downregulated_genes2 <- sum(df2$log2FoldChange < -1 & df2$padj < 0.05)

EnhancedVolcano(df2,
                lab = df2$Row.names,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1)

#### Heatmap Construction ####

top501 <- res1_annotated_filtered[1:50,1]
top501_expression <- assay(vsd[top501])
top501_expression  <- top501_expression - rowMeans(top501_expression)

res1_annotated_filtered_filter <- as.data.frame(res1_annotated_filtered)
res1_annotated_filtered_filter <- res1_annotated_filtered_filter %>% 
  filter(Row.names == row.names(top501_expression))

write.csv(
  res1_annotated_filtered_filter,
  "C:/Users/felip/Documents/Mestrado/Tese/Results/Transport/4VScontrol/Heatmap_genes_0vs4.csv",
  sep="\t", row.names = FALSE
)

h1 <- pheatmap(
  top501_expression,
  cluster_cols = FALSE,
  fontsize_row = 12,
  fontsize_col = 10
)

ggsave("../Results/Transport/HeatmapStress.tiff", h1,
       width = 26*1/2, height = 19.5*1/2, dpi = 300)

#### Heatmap for 24h Recovery ####

top502 <- res2_annotated_filtered[1:50,1]
top502_expression <- assay(vsd[top502])
top502_expression  <- top502_expression - rowMeans(top502_expression)

h2 <- pheatmap(
  top502_expression,
  fontsize_row = 12,
  cluster_cols = FALSE,
  fontsize_col = 10
)

ggsave(
  "C:/Users/felip/Documents/Mestrado/Tese/Results/Transport/24VScontrol/New_Analysis/HeatmapRecovery_none.tiff",
  h2, width = 26*1/2, height = 19.5*1/2, dpi = 300
)

#### Annotated Heatmaps with Condition Colors ####

annot <- as.data.frame(colData(vsd))
annot["Condition"] <- c(rep("Control_0h", 4),
                        rep("Stress_4h", 4),
                        rep("Control_24h", 4),
                        rep("Recovery_24h", 4))

anna_color = c(
  Control_0h = "darkgreen",
  Control_24h = "beige",
  Recovery_24h = "olivedrab",
  Stress_4h ="tan"
)

heat1 <- pheatmap(
  top501_expression,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annot["Condition"],
  annotation_colors = list(Condition = anna_color)
)

heat2 <- pheatmap(
  top502_expression,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annot["Condition"],
  annotation_colors = list(Condition = anna_color)
)

#### PlotCounts for Individual Gene ####

dds$Temperature <- as.factor(dds$Treatment)
dds$Time <- as.factor(dds$Time)

plotCounts(dds, gene= res1_annotated_filtered[1,1], intgroup=c("Treatment","Time"))

d <- plotCounts(dds, gene=res1_annotated_filtered[1,1],
                intgroup=c("Treatment","Time"), returnData = TRUE)

ggplot(d, aes(x=Treatment, y=count, col=Time)) +
  geom_point(position=position_jitter(w=0.1,h=0), size=1.5) +
  scale_y_log10() +
  theme_bw()
