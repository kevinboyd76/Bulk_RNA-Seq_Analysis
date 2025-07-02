# load libraries
library(DESeq2)
library(tidyverse)
library(vsn)
library(pheatmap)
library(clusterProfiler)
library(org.Dr.eg.db)
library(biomaRt)
library(ggrepel)

################################################################################
##                     Load the data and setup dds object                     ##
################################################################################
# set working directory
setwd("/Volumes/scratch/cobre-dev-bio/boydk/Varshney_Lab/DELETE_ME2/results/counts")

# read in the counts data, set rownames, and setup coldata object
counts_data = read.csv("merged_counts.csv")
rownames(counts_data) = counts_data[,1]
counts_data = counts_data[,-1]
sample_names = colnames(counts_data)
condition = factor(c("mutant", "mutant", "mutant", "control", "control", "control"))
colData = data.frame(row.names = sample_names, condition = condition)

# construct a DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ condition)
dds
no_filter = nrow(dds)
message("Genes retained before filtering: ", nrow(dds))

# filtering 
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]
dds
filter = nrow(dds)
message("Genes retained after filtering: ", nrow(dds))
message("Genes removed after filtering: ", no_filter - filter)
# factor level
dds$condition = relevel(dds$condition, ref = "control")

# Run DESeq2
dds = DESeq(dds)
res = results(dds)
res

################################################################################
##                             Observe the Results                            ##
################################################################################

# explore results
summary(res)

# change adjusted pvalue cutoff
res = results(dds, alpha = 0.05)
summary(res0.05)

# MA plot
plotMA(res)

# count plots
plotCounts(dds, gene = which.min(res$padj), intgroup = "condition")


################################################################################
# Extracting Transformed Values
vsd <- vst(dds, blind=FALSE)
# PCA plot
plotPCA(vsd, intgroup = "condition")
head(assay(vsd), 3)
meanSdPlot(assay(vsd))


################################################################################
##                               Volcano Plots                                ##
################################################################################

# Convert to data frame
res_df <- as.data.frame(res)

# Optional: filter out NA rows
#res_df <- res_df[!is.na(res_df$padj), ]

################################################################################
# Basic volcano plot

vol_plot = ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-log10 Adjusted P-value")
vol_plot


################################################################################

# Add a significance label
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
vol_plot2 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "gray")) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value",
       color = "Significance")
vol_plot2

################################################################################

#  Add a gene column from rownames
res_df$gene <- rownames(res_df)

# Remove rows with NA padj
res_df_filtered <- res_df[!is.na(res_df$padj), ]

# Pick top 30 by padj value
top_hits <- res_df_filtered %>%
  arrange(padj) %>%
  slice(1:40)

# 4. Make the volcano plot with labels
vol_plot <- ggplot(res_df_filtered, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4, aes(color = padj < 0.05 & abs(log2FoldChange) > 1)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_text_repel(data = top_hits, aes(label = gene), size = 4, max.overlaps = Inf) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = "Significant"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="none")

vol_plot


################################################################################
##                                 Heat Plots                                 ##
################################################################################
# Highest variance genes (not necessarily DEGs)

# Get top 30 most variable genes
top_genes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)

# Draw heatmap
pheatmap(assay(vsd)[top_genes, ],
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE, annotation_col = colData)

################################################################################
# All Genes Variance

# VST-transformed matrix
vsd_mat <- assay(vsd)

# Optional: scale across genes to center & normalize
vsd_scaled <- t(scale(t(vsd_mat)))

# Heatmap of all genes, no row labels
pheatmap(vsd_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = colData,
         fontsize_col = 10,
         main = "Heatmap of All Genes (VST-scaled)")

################################################################################
# Most statistically significant genes, even if fold change is small

# Remove rows with NA padj
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

# Get top 30 DEGs by adjusted p-value
top_degs <- rownames(res_df[order(res_df$padj), ])[1:30]

# Draw heatmap for those genes
pheatmap(assay(vsd)[top_degs, ],
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = colData)

################################################################################

# Statistically significant and biologically meaningful DEGs

# Filter for DEGs with padj < 0.05 and abs(log2FC) > 1
deg <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

deg2 = deg[,c(1,2,6)]
deg2_ordered <- deg2[order(deg2$log2FoldChange, decreasing = TRUE), ]


write.table(deg2_ordered, file = "Significant_DEGs.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

# Get top 30 from this filtered list
top_degs <- rownames(deg[order(deg$padj), ])[1:30]

# Draw heatmap for those genes (ordered by DE & Significance)
pheatmap(assay(vsd)[top_degs, ],
         cluster_rows = FALSE, cluster_cols = FALSE,
         show_rownames = TRUE,
         annotation_col = colData)

# Draw heatmap for those genes
pheatmap(assay(vsd)[top_degs, ],
         cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = TRUE,
         annotation_col = colData)
################################################################################
##                              Pathway Analysis                              ##
################################################################################

# pathway analysis
#res <- results(dds)
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$padj), ]

# Filter DEGs (adjust thresholds as needed)
deg <- res_df[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ]

# Connect to Ensembl for Danio rerio
mart <- useMart("ensembl", dataset = "drerio_gene_ensembl")

# Convert gene symbols to Entrez IDs
gene_symbols <- rownames(deg)
mapping <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                 filters = "external_gene_name",
                 values = gene_symbols,
                 mart = mart)

# Merge with DEG table
deg$gene <- rownames(deg)
deg_merged <- merge(deg, mapping, by.x = "gene", by.y = "external_gene_name")



# Create vector of Entrez IDs
entrez_genes <- deg_merged$entrezgene_id

# KEGG
kegg_res <- enrichKEGG(gene = entrez_genes,
                       organism = "dre",
                       pvalueCutoff = 0.05)

# GO Biological Process
go_res <- enrichGO(gene = entrez_genes,
                   OrgDb = org.Dr.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

# Visualize
barplot(go_res, showCategory = 20, title = "GO Biological Process")
dotplot(kegg_res, showCategory = 20, title = "KEGG Pathways")
