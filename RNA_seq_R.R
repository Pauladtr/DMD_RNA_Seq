# Load necessary packages
library(edgeR)
library(ggplot2)
library(pheatmap)
library (limma)
library (readxl)

# Set the work directory
setwd("C:/Users/paula/OneDrive/Documents/EdgeR_DMD")

# Read in raw count data
counts <- read.table("countData.txt", header = TRUE, row.names = 1)

# Read sample information
sampleInfo <- read_excel("design.csv.xlsx", sheet = 1)

#Creating the DGE object
dgeFull <- DGEList(counts, group=sampleInfo$condition)
dgeFull


# Filter out lowly expressed genes

#Transformations
cpm <- cpm(dgeFull)
lcpm <- cpm(dgeFull, log = TRUE)

#Remove genes that are lowly expressed
keep <- filterByExpr(dgeFull, group = sampleInfo$Condition, min.count = 10)
dgeFilt <- dgeFull[keep, ]
D_Filtered <- dgeFull[keep.exprs, ]

# Normalization of gene expression
D_Normalized <- calcNormFactors(D_Filtered)


# Create a DGEList object
DGE_List <- DGEList(counts = D_Normalized, group = D_Filtered$samples$group)


# Create MDS plot

# Create a vector of colors based on sample names
colors <- c(rep("red", 3), rep("blue", 3), rep("green", 3))

# Assign colors to the DGE_List object
DGE_List$samples$color <- colors

# Create MDS plot with colored groups
mds <- plotMDS(DGE_List, col = DGE_List$samples$color, main = "MDS plot")

