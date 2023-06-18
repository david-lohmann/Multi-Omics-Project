library(tximport)

# First, set the R working directory to the folder containing Kallisto results 
setwd("") #Complete with the path in your computer

# Load the file containing the EnsemblIDs and gene names
load('gencode.v38_geneInfo.RData');
head(geneInfo)

# Define filenames to import
folders <- dir(pattern="rep")
fileNames <- file.path(".", folders, 'abundance.tsv');
names(fileNames) <- folders 

# Import all data using tximport function
readCounts_gene <- tximport(fileNames, type = 'kallisto',txOut = F,tx2gene=geneInfo[,1:2], ignoreTxVersion = T,countsFromAbundance = "lengthScaledTPM")

# Explore the object created with tximport
names(readCounts_gene)
View(head(readCounts_gene$abundance)) # matrix containing TPMs
View(head(readCounts_gene$counts)) # matrix containing read counts
head(readCounts_gene$length) 
head(readCounts_gene$countsFromAbundance) # indicates the type of normalization applied to get the values in "abundance"

# To facilitate interpretation of downstream results, replace EnsemblIDs with Official gene symbol (i.s gene name)
geneNames <- geneInfo[match(rownames(readCounts_gene$abundance), geneInfo$Ensembl_GeneID),"GeneSymbol"]
rownames(readCounts_gene$abundance) <- rownames(readCounts_gene$counts) <- geneNames
head(readCounts_gene$abundance)
head(readCounts_gene$counts)

# Save readCounts in an RData file for downstream analyses
save(readCounts_gene, file="readCounts_gene.RData")


library(PCAtools)
# Prepare gene expression matrix for PCA analysis: Step 1) get log TPMs (we add 1 TPM to each gene to avoid infinite values after log)
logTPMs <- log2(readCounts_gene$abundance+1)

# Prepare gene expression matrix for PCA analysis: Step 2) remove duplicated genes
uniqueGenes <- unique(rownames(logTPMs))
logTPMs <- logTPMs[uniqueGenes,]

#Prepare metadata with sample type
sampleTypes <- gsub("_rep[1234]", "", colnames(logTPMs))
metaData <- data.frame(sampleTypes); rownames(metaData) <- colnames(logTPMs)
metaData

#Run PCA
pca.res <- pca(logTPMs, metadata=metaData)

#Plot variance explained by each component
screeplot(pca.res)

#Plot 2 selected components/eigenvectors
biplot(pca.res)
biplot(pca.res, colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with colors by sample type
biplot(pca.res, lab="", colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot without sample names
biplot(pca.res, x="PC2", y="PC3",lab="",colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top') # Biplot with PC1 and PC3

#Plot several components
pairsplot(pca.res, colby="sampleTypes")
pairsplot(pca.res)

# Plot the component loadings and label genes most responsible for variation
PC_genes <- pca.res$loadings
PC1_genes <- PC_genes[order(PC_genes$PC1, decreasing=T),]
head(PC1_genes)
tail(PC1_genes)

plotloadings(pca.res, components = c("PC1", "PC2", "PC3"),rangeRetain =0.1) #retaining 1% of the loadings per PC

# Produce barplot to confirm expression levels of genes associated with Principal Component 1
#plotCol <- rep(c("steelblue", "gray"), each=3)
plotCol <- ifelse(grepl("wt",sampleTypes),"steelblue","gray")
barplot(logTPMs["MYBL2",], col=plotCol, las=2, main="MYBL2", ylab="Expression levels (logTPMs)") #Gene positively correlated with PC1
legend("topleft", fill=unique(plotCol), legend=c("MYBL2","ctrl"), bty="n")
barplot(logTPMs["ATP6V0E2",], col=plotCol, las=2, main="SNORD89", ylab="Expression levels (logTPMs)") #Gene negatively correlated with PC1
legend("topright", fill=unique(plotCol), legend=c("actFOXO1","ctrl"), bty="n")
#SNORD89,MIR4292,RNU6-216P
#RDH11, ATP6V0E2


pdf("PCA_plots.pdf") # open pdf file
biplot(pca.res, colby="sampleTypes", hline = 0, vline = 0,legendPosition = 'top')
plotloadings(pca.res, components = c("PC1", "PC2", "PC3"),rangeRetain =0.1) 
dev.off()



library(edgeR)

# Create DGEList data class (specific for edgeR package)
y <- DGEList(counts=readCounts_gene$counts, group= sampleTypes)


# Filter out lowly expressed genes
keep <- filterByExpr(y, group=sampleTypes)
y <- y[keep, ,keep.lib.sizes=FALSE]

# Normalize for library sizes
y <- calcNormFactors(y)
y$samples


# Define the design and contrast matrix based on the experimental design, meaning define which comparison to be made
design_matrix <-ifelse(grepl("wt",sampleTypes),"wt","knockout")
design_matrix <- model.matrix(~0+sampleTypes)
colnames(design_matrix) <- gsub("sampleTypes", "", colnames(design_matrix))
rownames(design_matrix) <- colnames(logTPMs)
design_matrix
#contrast_matrix <- makeContrasts(wt-ctrl, levels = design_matrix)
contrast_matrix <- makeContrasts(kdMYBL2-WT, levels = design_matrix)
contrast_matrix


# Estimate the dispersion (Biological coefficient of variation) and Fit model
y <- estimateDisp(y, design_matrix)
plotBCV(y)
fit <- glmQLFit(y,design_matrix) # Fit a quasi-likelihood negative binomial generalized log-linear model to count data. 


# Differential expression analysis
qlf <- glmQLFTest(fit, contrast =  contrast_matrix) # test for differential expression between the experimental groups using quasi-likelihood F-test

# Get differentially expressed genes (DEGs)
summary(decideTests(qlf,p.value = 0.05,adjust.method = "fdr"))
DEGs <- topTags(qlf, nrow(qlf$table), p.value=0.05, adjust.method = "fdr")$table
upReg <- rownames(DEGs)[which(DEGs$logFC > 0)]
downReg <- rownames(DEGs)[which(DEGs$logFC < 0)]
length(upReg)
length(downReg)
# Write table with DEGs
library(openxlsx)
write.xlsx(DEGs, file = "DEGs.xlsx", rowNames=T)

# Produce Volcano Plot (customized)
allGenes <- topTags(qlf, n = nrow(qlf$table), p.value = 1)$table
plotData <- cbind(allGenes$logFC, -log10(allGenes$FDR)); rownames(plotData) <- rownames(allGenes)
plot(plotData, pch=20,col="gray",xlab="Biological Variation (log2 Fold-Change)", ylab="Statistical Significance (-log10 P-Value)")
abline(h=-log10(0.05), v=c(-1,1),lty=2)
points(plotData[upReg,], col="tomato", pch=20)
points(plotData[downReg,], col="steelblue", pch=20)
text(plotData[upReg[1:5],], labels=upReg[1:5],col="tomato", pos=sample(c(1:3), size=10, replace=T), cex=0.8)
text(plotData[downReg[1:5],], labels=downReg[1:5],col="steelblue", pos=sample(c(1:2), size=10, replace=T), cex=0.8)

# Produce MA plot (customized)
plotData <- cbind(allGenes$logCPM, allGenes$logFC); rownames(plotData) <- rownames(allGenes)
plot(plotData, pch=20,col="gray",xlab="Mean Expression Levels (log CPMs)", ylab="Biological Variation (log2 Fold-Change)")
abline(h=c(-1,1),lty=2)
points(plotData[upReg,], col="tomato", pch=20)
points(plotData[downReg,], col="steelblue", pch=20)
text(plotData[upReg[1:5],], labels=upReg[1:5],col="tomato", pos=sample(c(1:3), size=10, replace=T), cex=0.8)
text(plotData[downReg[1:5],], labels=downReg[1:5],col="steelblue", pos=sample(c(1:2), size=10, replace=T), cex=0.8)


# Heatmap with Top DEGs
library(gplots); library(RColorBrewer)
plotCol_exp <- brewer.pal(9, "Purples")
heatmap.2(logTPMs[c(upReg[1:10], downReg[1:10]),], scale="row", trace="none", density.info="none", ColSideColors = plotCol, col=plotCol_exp)


# Get all expressed genes (DEGs)
summary(decideTests(qlf,p.value = 0.05,adjust.method = "fdr"))
EXPGenes <- topTags(qlf, nrow(qlf$table))$table
upReg <- rownames(DEGs)[which(DEGs$logFC > 0)]
downReg <- rownames(DEGs)[which(DEGs$logFC < 0)]

# Write table with DEGs
library(openxlsx)
#write(EXPGenes, file = "EXPGenes.rnk")
write.xlsx(EXPGenes, file = "EXPGenes.xlsx", rowNames=T)


EXPGenes$rank <- (sign(EXPGenes$logFC)* - log10(EXPGenes$PValue)) 
head(EXPGenes)

names_ranks = data.frame(rownames(EXPGenes),EXPGenes$rank)
head(names_ranks)



write.table(names_ranks, file = "expressed_genes_ranks.rnk", append = FALSE, sep = "\t", row.names = F, col.names = FALSE, quote = F,dec = ".")
