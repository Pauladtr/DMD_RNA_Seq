# Load necessary packages
library(edgeR)
library(ggplot2)
library(pheatmap)
library (limma)
library (readxl)
library(RColorBrewer)
library(gplots)

# Set the work directory
setwd("C:/Users/paula/OneDrive/Documents/EdgeR_DMD")

# Read in raw count data
counts <- read.table("countData.txt", header = TRUE, row.names = 1)

# Read sample information
sampleInfo <- read_excel("design.csv.xlsx", sheet = 1)

#Creating the DGE object
dgeFull <- DGEList(counts, group=sampleInfo$condition)
dgeFull

#Transformations
cpm <- cpm(dgeFull)
lcpm <- cpm(dgeFull, log = TRUE)

##Create a density plot using ggplot2 (before filtration)

# Convert lcpm to a numeric vector
lcpm <- as.vector(lcpm)
ggplot(data.frame(lcpm), aes(x = lcpm)) +
  geom_density() +
  labs(x = "lcpm", y = "Density") +
  ggtitle("Density Plot")

##Filter out lowly expressed genes

#Remove genes that are lowly expressed
keep <- filterByExpr(dgeFull, group = sampleInfo$Condition, min.count = 10)
dgeFilt <- dgeFull[keep, ]
D_Filtered <- dgeFull[keep, ]


## Create a density plot using ggplot2 (after filtration)

# Convert lcpm to a numeric vector
lcpm_filtered <- as.vector(cpm(dgeFilt, log = TRUE))
ggplot(data.frame(lcpm_filtered), aes(x = lcpm_filtered)) +
  geom_density() +
  labs(x = "lcpm", y = "Density") +
  ggtitle("Density Plot after Filtration")

#Transformations
cpm <- cpm(dgeFull)
lcpm <- cpm(dgeFull, log = TRUE)

#Density plot - Raw Data
nsamples <- ncol(dgeFull)

#Make a color vector
colors <- c(rep("red", 3), rep("blue", 3), rep("green", 3))

par(mfrow = c(1, 2))
plot(density(lcpm[, 1]), ylim = c(0, 1.35), col = colors[1], lwd = 2, las = 2, 
     main = "Raw Data", xlab = "logCPM")
for(i in 2:nsamples){
  den <- density(lcpm[, i])
  lines(den$x, den$y, lwd = 2, col = colors[i])
}

#Remove low expressed genes
keep.exprs <- filterByExpr(dgeFull, group = sampleInfo$condition, min.count = 10)
dgeFiltered <- dgeFull[keep.exprs, , keep.lib.sizes = FALSE]

#Transformations
cpm <- cpm(dgeFiltered)
lcpm <- cpm(dgeFiltered, log = TRUE)

#Density plot - Filtered Data
nsamples <- ncol(dgeFiltered)
plot(density(lcpm[, 1]), ylim = c(0, 1.35), col = colors[1], lwd = 2, las = 2, 
     main = "Filtered Data", xlab = "logCPM")
for(i in 2:nsamples){
  den <- density(lcpm[, i])
  lines(den$x, den$y, lwd = 2, col = colors[i])
}
# Normalization of gene expression
D_Normalized <- calcNormFactors(D_Filtered)


## Create MDS plot

# Create a vector of colors based on sample names
colors <- c(rep("red", 3), rep("blue", 3), rep("green", 3))

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
logcounts <- cpm(dgeFilt, log=TRUE)

# Boxplot after normalization
boxplot(logCounts, 
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

