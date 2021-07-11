library(EnhancedVolcano)
library(edgeR)
library(statmod)
library(gmp)

setwd('Your Directory')

######Calculate Differentially Expressed Genes in Non-SPG Cell Types Between Diff vs Undiff SPG Neighborhood##############

# Read in the neighborhood gene expression matrix outputted from the SPG_Compartment_Analysis.ipynb. 
# Note that in order to increase statistical power, we combined data from three replicates (i.e., 3 Slide-seq arrays).
counts <- read.csv("Combined_DGE_SPG_Neighboring_Beads_K5.csv",row.names=1)
counts[1:5, 1:3]

#Filter lowly-expressed genes
countdata <- counts[ rowSums(counts) >= 90, ]
dim(countdata)

# Read in the combined meta-data file
coldata <- read.csv("Combined_Metadata_SPG_Neighboring_Beads_K5.csv",row.names=1)
coldata[1:5, 1:4]

#Change the data type of the Sub_cell_type column from integer to character
coldata[,'Sub_cell_type']<-factor(coldata[,'Sub_cell_type'])
type(coldata['Sub_cell_type'])

#Check if the sample names in the meata data match the sample names in the count data
all(rownames(coldata) %in% colnames(countdata))

#Check if the samples are arranged in the same order in meta data and the count data
all(rownames(coldata) == colnames(countdata))

#Reorder the samples in the count data if not the same
countdata <- countdata[, rownames(coldata)]

#Check if the samples are arranged in the same order in meta data and the count data
all(rownames(coldata) == colnames(countdata))

#Create an edgeR object
allcounts = DGEList(counts=countdata)

#Calculate normalization factors. All samples should have similar values
allcounts <- calcNormFactors(allcounts)
allcounts$samples

#Design matrix for edgeR. The last variable is always the one of interest
design <- model.matrix(~ Replicate+Sub_cell_type+Neighborhood, data=coldata)
rownames(design) <- colnames(allcounts)
head(design)

allcounts = estimateDisp(allcounts, design, robust=TRUE)

#This shows the curve fitting to reestimate the dispersion
plotBCV(allcounts)

#The recommended fitting procedure in edgeR returns Quasi-F as statistic for significance in the GLM
fit = glmQLFit(allcounts, design, robust=TRUE)
lrt = glmQLFTest(fit)
plotMD(lrt)

# The topTags function gives the genes where the variable was significant
# if we want all genes, we need to ask for it using the n parameter
topgenes<-topTags(lrt,n=Inf)

# Check if there is any differentially expressed gene 
is.de <- decideTestsDGE(lrt)
summary(is.de)

res = topgenes@.Data[[1]]
summary(res)

write.csv(res, "DGE_Results.csv")
