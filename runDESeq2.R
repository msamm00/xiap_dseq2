#load libraries

library(DESeq2)
library(tidyverse)
library(airway)
library(ggplot2)
library(pcaExplorer)
library(pheatmap)
library(AnnotationDbi)
library("org.Hs.eg.db")
library("EnsDb.Hsapiens.v86")
library("ggfortify")
citation("pcaExplorer")
options(ggrepel.max.overlaps = Inf)
# Step 1: preparing count data ----------------

# read in counts data
counts_data <-read.csv("./GSE151648_liver-iri-counts_cleaned.txt", sep="", stringsAsFactors=TRUE)
# read in sample info
colData <- read.csv("./SraRunTable.txt")
# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))
# are they in the same order?
all(colnames(counts_data) == rownames(colData))

#Mapping geneids to gene symbols /names for all genes
symbol<- mapIds(EnsDb.Hsapiens.v86,keys= rownames(counts_data),keytype="GENEID",column="SYMBOL", multiVals="first")

for (row in 1:nrow(counts_data))  {
  tryCatch (
    expr = {      
    rownom <- rownames(counts_data[row,])
    name <- ifelse(is.na(symbol[rownom]), rownom, symbol[rownom])
    row.names(counts_data)[row] <- name},
    error = function(e) {  
     },
    finally = {  
    }
  )
}

# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = (counts_data),
                              colData = colData,
                              design = ~ iri+Timepoint)
dds
# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# set the factor level
dds$iri <- relevel(dds$iri, ref = "IRIMinus")

# NOTE: collapse technical replicates if present
# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)
summary(res)
#Regularized log transformation
vst <- varianceStabilizingTransformation(dds,blind=TRUE,fitType = "parametric" ) 

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(vst)


#Understand differential expressed genes using pcaexplorer 
gene_annotations <- data.frame(row.names=rownames(counts_data),"gene_name"=ifelse(is.na(symbol), rownames(counts_data), symbol))
pcaExplorer(dds = dds, dst=vst,annotation=gene_annotations)


#boxplot genes of interest
normalizedCounts <- counts(dds, normalized=T)


#PCA plot on data sample 
z<-plotPCA(vst, intgroup=c("iri","Timepoint"))
z

#Get counts data for a given gene of interest
Counts <- plotCounts(dds, gene="XIAP", intgroup=c("iri","Timepoint"), normalized=TRUE,returnData=TRUE)

#PCA box plot data across iri and timepoint 
ggplot(Counts, aes(x=iri, y=count,colour=Timepoint)) + 
  geom_jitter(width = 0.1) +
  geom_point(position=position_jitter(w=0.1,h=0.1)) + 
  scale_y_log10() + geom_boxplot()
#Following stimulation with LPS, Kupffer cells produce numerous inflammatory cytokines,
#such as TNF-α, IL-1β, IL-6, IL-12, IL-18, IL-10, and several chemokines [42]. Kupffer cell-derived IL-12 and IL-18 
#activate hepatic natural killer (NK) cells to increase the synthesis and release of antimicrobial IFN-γ [43].

#heatmaps to compare correlation across genes
#rowsOfInterest <- c("BIRC2","BIRC3","TNFAIP3","XIAP","BIRC5","BIRC6","BIRC7" ,"IL6","CXCL8","RIPK1","DIABLO","BCL2","IL1B","TNF","TLR4")
rowsOfInterest <- c("BCL2","TNF","CXCL8","IL6")
mat  <- assay(vst)[rowsOfInterest,]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vst)[, c("iri","Timepoint")])
#Column/global scaling to compare across genes
pheatmap(mat,annotation=anno,scale="column")



