# Load necessary packages
library(edgeR)
library(ggplot2)
library(pheatmap)
library (limma)
library (readxl)
library(RColorBrewer)
library(gplots)
library (reshape2)
library(magrittr)
library(dplyr)

# Set the work directory
setwd("C:/Users/paula/OneDrive/Documents/EdgeR_DMD")

# Read in raw count data
counts <- read.table("countData.txt", header = TRUE, row.names = 1)

# Read sample information
sampleInfo <- read_excel("design.csv.xlsx", sheet = 1)

#Creating the DGE object
dgeFull <- DGEList(counts, group=sampleInfo$Condition)
dgeFull

#Transformations
lcpm <- cpm(dgeFull, log = TRUE)

##Create a density plot using ggplot2 (before filtration)
# Create a vector of colors based on sample names
colors <- c(DMDII = "red", DMDI = "blue", WT = "green")

ggplot(lcpm %>% reshape2::melt() %>%
         mutate(group = gsub("_.*","",Var2)), aes(x = value, group = group, colour = group)) +
  geom_density() +
  theme_bw() +
  labs(x = "lcpm", y = "Density") +
  scale_color_manual(values = colors) +
  ggtitle("Density Plot")

##Filter out lowly expressed genes

#Remove genes that are lowly expressed
keep <- filterByExpr(dgeFull, group = sampleInfo$Condition, min.count = 40)
D_Filtered <- dgeFull[keep, ]


## Create a density plot using ggplot2 (after filtration)
ggplot(reshape2::melt(cpm(D_Filtered, log = TRUE)) %>% 
         mutate(group = gsub("_.*", "", Var2)), 
       aes(x = value, group = group, colour = group)) +
  geom_density() +
  theme_bw() +
  labs(x = "lcpm", y = "Density") +
  scale_color_manual(values = colors) +
  ggtitle("Density Plot (After Filtration)")


# Normalization of gene expression
D_Normalized <- calcNormFactors(D_Filtered)


## Create MDS plot

# Assign colors to the DGE_List object
dgeFilt$samples$color <- colors

# Create MDS plot with colored groups
mds <- plotMDS(dgeFilt, col = dgeFilt$samples$color, main = "MDS plot")


## Boxplot
# Extract the counts from dgeFull
counts_before_norm <- dgeFull$counts

# Boxplot before the normalization
boxplot(counts_before_norm, 
        main = "Boxplot before Normalization", 
        xlab = "", 
        ylab = "Log2 counts per million",
        las = 2,
        col = D_Filtered$samples$color,
        cex.axis = 0.7)

# Check distributions of samples using boxplots
lpcm <- cpm(dgeFilt, log=TRUE)

# Boxplot after normalization
boxplot(lpcm, 
        xlab="", 
        ylab="Log2 counts per million",
        las=2,
        col = D_Filtered$samples$color,  # Use the color vector from the MDS plot
        cex.axis = 0.7,  # Adjust the size of the sample names
        mar = c(5, 6, 4, 2),  # Adjust the plot margins 
        main= "Boxplot after Normalization"
)

##Heatmap

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(lcpm, 1, var)
head(var_genes)

# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
head(select_var)

# Subset logcounts matrix
highly_variable_lcpm <- lcpm[select_var,]
dim(highly_variable_lcpm)

head(highly_variable_lcpm)

# Get some colours
mypalette <- brewer.pal(11, "RdYlBu")

# http://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3
morecols <- colorRampPalette(mypalette)

# Set up colour vector 
col.cell <- c("purple", "orange")[sampleInfo$Condition]

# Plot the heatmap
dev.new()
heatmap.2(filtered_normalized_selected,
          col = rev(morecols(50)),
          trace = "column",
          main = "Top 500 most variable genes across samples",
          ColSideColors = col.cell,
          cexCol = 0.9,
          scale = "row", 
          labRow = FALSE,
          )

##Differential Gene expression analysis
dgeFull <- estimateCommonDisp(D_Normalized)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull

#Define the groups for comparison
group1 <- "DMDI"
group2 <- "DMDII"
group3 <- "WT"

#Perform an exact test for the difference between the conditions: DMDI and WT, DMDII and WT
dgeTest1 <- exactTest(dgeFull, pair = c(group1, group3))
dgeTest1
dgeTest2 <- exactTest(dgeFull, pair = c(group2, group3))
dgeTest2

# Get the differential expression results for each comparison
de_genes1 <- topTags(dgeTest1, n = Inf)$table
de_genes2 <- topTags(dgeTest2, n = Inf)$table

# Print the differential expression results
head(de_genes1)
head(de_genes2)

#Get a summary DGE table (returns significant genes with absolute log fold change at least 1 and adjusted p value < 0.05)
summary(decideTests(object = dgeTest1, lfc = 1))
summary(decideTests(object = dgeTest2, lfc = 1))

#Export differential gene expression analysis table to CSV file
write.csv(as.data.frame(de_genes1), file="condition_DMDI_vs_WT_dge.csv")
write.csv(as.data.frame(de_genes2), file="condition_DMDII_vs_WT_dge.csv")

# Merge DMDI and DMDII results into one table
merged_results <- merge(de_genes1, de_genes2, by = "row.names", all = TRUE)
colnames(merged_results) <- c("Gene", "DMDI_logFC","DMDI_logPCM", "DMDI_PValue", "DMDI_FDR", "DMDII_logFC","DMDII_logPCM", "DMDII_PValue", "DMDII_FDR")
merged_results <- merged_results[order(merged_results$Gene), ]  # Sort by gene name

# Export merged results to CSV file
write.csv(merged_results, file = "merged_results.csv", row.names = FALSE)

##DGE vizualisation
# Volcano plot for de_genes1 (DMDI vs. WT)
ggplot(de_genes1, aes(x = logFC, y = -log10(PValue))) +
  geom_point(size = 1, color = "black") +
  geom_point(data = subset(de_genes1, abs(logFC) >= 1 & PValue < 0.05), size = 1, color = "red") +
  theme_bw() +
  labs(x = "Log Fold Change (DMDI vs. WT)", y = "-log10(P-value)") +
  ggtitle("Volcano Plot (DMDI vs. WT)")

# Volcano plot for de_genes2 (DMDII vs. WT)
ggplot(de_genes2, aes(x = logFC, y = -log10(PValue))) +
  geom_point(size = 1, color = "black") +
  geom_point(data = subset(de_genes2, abs(logFC) >= 1 & PValue < 0.05), size = 1, color = "red") +
  theme_bw() +
  labs(x = "Log Fold Change (DMDII vs. WT)", y = "-log10(P-value)") +
  ggtitle("Volcano Plot (DMDII vs. WT)")



#MA plot for dgeTest1 and degTest2

plotMD(dgeTest1)
plotMD(dgeTest2)

##Finding the overlap between the de_genes1 and de_genes2

# Filter de_genes1 for upregulated genes
upregulated_genes1 <- de_genes1[de_genes1$logFC > 1 & de_genes1$FDR < 0.05, ]

# Filter de_genes2 for upregulated genes
upregulated_genes2 <- de_genes2[de_genes2$logFC > 1 & de_genes2$FDR < 0.05, ]

# Get the gene names for upregulated genes
upregulated_gene_names1 <- rownames(upregulated_genes1)
upregulated_gene_names2 <- rownames(upregulated_genes2)

# Find the intersect of upregulated genes
intersect_upregulated_genes <- intersect(upregulated_gene_names1, upregulated_gene_names2)

# Print the intersecting upregulated genes
print(intersect_upregulated_genes)



## Filter de_genes1 for downregulated genes
downregulated_genes1 <- de_genes1[de_genes1$logFC < -1 & de_genes1$FDR < 0.05, ]

# Filter de_genes2 for downregulated genes
downregulated_genes2 <- de_genes2[de_genes2$logFC < -1 & de_genes2$FDR < 0.05, ]

# Get the gene names for downregulated genes
downregulated_gene_names1 <- rownames(downregulated_genes1)
downregulated_gene_names2 <- rownames(downregulated_genes2)

# Find the intersect of downregulated genes
intersect_downregulated_genes <- intersect(downregulated_gene_names1, downregulated_gene_names2)

# Print the intersecting downregulated genes
print(intersect_downregulated_genes)

# Find genes upregulated in de_genes1 but downregulated in de_genes2
genes_upregulated_de1_downregulated_de2 <- intersect(upregulated_gene_names1, downregulated_gene_names2)
# Print the list of genes
print(genes_upregulated_de1_downregulated_de2)

#Find genes downregulated in de_genes1 but upregulated in de_genes2
genes_downregulated_de1_upregulated_de2 <- intersect(downregulated_gene_names1, upregulated_gene_names2)
# Print the list of genes
print(genes_downregulated_de1_upregulated_de2)

##Save the lists 

# Export the intersecting genes as Excel files
write.csv(as.data.frame(intersect_upregulated_genes),
           file = "intersect_upregulated_genes.csv")

write.csv(as.data.frame(intersect_downregulated_genes),
           file = "intersect_downregulated_genes.csv")

write.csv(as.data.frame(genes_upregulated_de1_downregulated_de2),
           file = "genes_upregulated_de1_downregulated_de2.csv")

write.csv(as.data.frame(genes_downregulated_de1_upregulated_de2),
           file = "genes_downregulated_de1_upregulated_de2.csv")

