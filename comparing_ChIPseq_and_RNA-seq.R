
peakfile <-read.table("MYBL2_peaks.narrowPeak")
colnames(peakfile) <- c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak")
head(peakfile)

# get a 400bp window around summit (-200 and +200)
Motif_Peaks <- data.frame((peakfile$chrom), (peakfile$chromStart-200), (peakfile$chromStart+200), stringsAsFactors = FALSE)
head(Motif_Peaks)

# create bed file
options(scipen=999)

write.table(Motif_Peaks, file= "Peaks_for_motif_detection.bed", 
            row.names=FALSE,col.names = FALSE, sep="\t", quote =FALSE)

# list files in directory and check if the file you created is there.
dir()



library("ChIPpeakAnno")
library("GenomicRanges")

#library("rtracklayer")
#bed <- import("/home/participant/Course_Materials/ChIPSeq/Materials/motifs/Peaks_for_motif_detection.bed")
#gr1 <- GRanges(seqnames=bed$V1, ranges=IRanges(start=bed$V2, end=bed$V3)) 

#read in the bed file
df1<-read.table("Peaks_for_motif_detection.bed", header=FALSE)

#convert the peaks to a GRanges object
gr1 <- GRanges(seqnames=df1$V1, ranges=IRanges(start=df1$V2, end=df1$V3))



##annotation
library("EnsDb.Hsapiens.v86")


## create annotation file from EnsDb (Ensembl) or TxDb (transcript annotation) packages

annoData <- toGRanges(EnsDb.Hsapiens.v86, feature="gene")
annoData[1:2]

# annotate the peak GRanges object with peaks mapped to gene with a -4000 and 500 bp window around the TSS

anno.gr1 <- annotatePeakInBatch(gr1, 
                                AnnotationData=annoData, 
                                output="nearestBiDirectionalPromoters",
                                bindingRegion=c(-4000, 500))

#trim out of bound ranges
anno.gr1 <- trim(anno.gr1)




#annotate with Entrez IDs

library("org.Hs.eg.db")
anno.gr1 <- addGeneIDs(anno.gr1,"org.Hs.eg.db",IDs2Add = "entrez_id")

# list annotated peaks
head(anno.gr1)

length(anno.gr1$gene_name)
#Número peaks


setwd("C:\\Users\\jdpl2\\OneDrive\\Ambiente de Trabalho\\Mestrado\\2º Semestre\\Métodos Computacionais em Multi-Ómicas\\Projeto\\data\\")
library(readxl)
DGEs <- read_excel("DEGs.xlsx")

commom_genes <- intersect(DGEs$...1,anno.gr1$gene_name)
length(commom_genes)
#Genes comum entre os peaks e os nossos DEGs

#No paper temos 10 genes que são regulado pelo MYBL2
sum(anno.gr1$gene_name == "CENPA") #0
sum(anno.gr1$gene_name == "FAM83D") #1
sum(anno.gr1$gene_name == "KNSTRN") #0
sum(anno.gr1$gene_name == "NEURL1B") #1
sum(anno.gr1$gene_name == "CCNF") #0
sum(anno.gr1$gene_name == "CDKN3") #1
sum(anno.gr1$gene_name == "CCNB1") #1
sum(anno.gr1$gene_name == "BORA") #1
sum(anno.gr1$gene_name == "CDK1") #1
sum(anno.gr1$gene_name == "CCNA2") #0
#Encontramos 6 genes que estão citados no paper

genes <- c("CENPA","FAM83D","KNSTRN","NEURL1B","CCNF","CDKN3","CCNB1","BORA","CDK1","CCNA2")

res <- DGEs[DGEs$...1 %in% genes, ] #Os genes estão todos presentes no DGEs

library(ggVennDiagram)

ggVennDiagram(list(DEGs = DGEs$...1,
                   'ChIP-seq genes' = anno.gr1$gene_name))
