# Referring to the source code https://github.com/MarioniLab/MNN2017/blob/master/Haematopoiesis/prepareData.R

# This script prepares data for the hematopoiesis analysis.
# It involves two publicly available datasets.

##########################################
##########################################
rm(list=ls())

# Setting the working directory
# setwd("D://Data_Set/BUSseq/BUSseq_implementation-master/Hematopoietic")

##############
## GSE81682 ## #<<<<<<<-------This is newly added.
##############
# Downloading and reading the counts, metadata of Nestorowa et al. 2016

# Downloading from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81682
fname <- "./RawData/GSE81682_HTSeq_counts.txt.gz"
if (!file.exists(fname)) { download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81682/suppl/GSE81682_HTSeq_counts.txt.gz", fname) }
dataF <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
dataF <- as.matrix(dataF)
dim(dataF)

fname <- "./RawData/metaF.txt"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/all_cell_types.txt", fname) }
metaF <- read.table(fname, stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
metainds <- match(colnames(dataF), rownames(metaF))
missing.meta <- is.na(metainds)
metaF <- metaF[metainds,] # This will contain NA's... which is okay, at this point, to preserve length.

# Defining the cell type based on the metadata.
metatypeF <- rep("other", nrow(metaF))
for (col in rev(colnames(metaF))) { # reverse, so earlier columns end up overwriting later ones.
  chosen <- metaF[,col]==1
  metatypeF[chosen] <- sub("[0-9]?_.*", "", col)
}
metatypeF[metatypeF=="ESLAM"] <- "HSPC"

# Filling in metadata from the cell sorting label, if metadata was missing.
metatypeF[missing.meta] <- sub("_.*", "", colnames(dataF)[missing.meta])
metatypeF[metatypeF=="LT-HSC"] <- "LTHSC"
metatypeF[metatypeF=="Prog"] <- "other"
colnames(dataF)<-metatypeF

# Cleaning up memory.
gc() 

##########################################
##########################################

##############
## GSE72857 ## #<<<<<<<-------This is newly added.
##############
# Download and read the counts and metadata of Paul et al. 2015

# Download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72857
fname <- "./RawData/GSE72857_umitab.txt.gz"
if (!file.exists(fname)) { download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72857/suppl/GSE72857_umitab.txt.gz", fname) }
dataA <- read.table(fname, header=TRUE, row.names=1)
metaA <- read.csv2("./RawData/MAP.csv",sep=",",stringsAsFactors = FALSE, head=TRUE, row.names=1)
dim(dataA)

# Only selecting cells that are in the metadata.
metainds <- match(rownames(metaA), colnames(dataA))
dataA <- dataA[,metainds]
dataA <- as.matrix(dataA)

# Organizing cell type labels.
metatypeA <- character(nrow(metaA))
metatypeA[metaA[,1]<7] <- "ERY"
metatypeA[metaA[,1]>6 & metaA[,1]<12] <- "CMP"
metatypeA[metaA[,1]>11] <- "GMP"
colnames(dataA) <- metatypeA

# Cleaning up memory.
gc()

##########################################
##########################################

# Downloading list of highly variable genes identified by Nestrowa et al. 2016
fname <- "./RawData/coordinates_gene_counts_flow_cytometry.txt.gz"
if (!file.exists(fname)) { download.file("http://blood.stemcells.cam.ac.uk/data/coordinates_gene_counts_flow_cytometry.txt.gz", fname) }
TFs <- read.table(fname, nrows=1, stringsAsFactors=FALSE)
features <- as.character(unlist(TFs))
features <- features[grep("ENSMUS", features)]

# Pulling down IDs from BioMaRt.
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="www.ensembl.org" )
out <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), values = features, mart = mart,filters = "ensembl_gene_id")

# Selecting features that are HVGs and present in both data sets.
mF <- match(out$ensembl_gene_id, rownames(dataF))
mA <- pmatch(out$mgi_symbol, rownames(dataA)) # partial, due to use of concatenated gene symbols.
# Three genes have been renamed in the BioMart database, so they fail to be matched my pmatch.
# We add them mutually
# "Dscr3" --> "Vps26c"
# referring to https://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000022898;r=16:94497783-94526830
print(out[112,])
mA[112] <- which(rownames(dataA)=="Dscr3")

# "Ssfa2" --> "Itprid2"
# referring to http://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000027007;r=2:79635352-79672966
print(out[1572,])
mA[1572] <- which(rownames(dataA)=="Ssfa2")

# "Fbxo18" --> "Fbh1"
# referring to http://asia.ensembl.org/Mus_musculus/Gene/Summary?g=ENSMUSG00000058594;r=2:11742573-11777582
print(out[4218,])
mA[4218] <- which(rownames(dataA)=="Fbxo18")

keep <- !is.na(mF) & !is.na(mA)

dataA2 <- dataA[mA[keep],]
dataF2 <- dataF[mF[keep],]
rownames(dataA2) <- rownames(dataF2)

if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}

save(dataA2,dataF2,file="./RawCountData/hemat_countdata.RData")

###########
# END