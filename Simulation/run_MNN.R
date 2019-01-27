#Apply MNN to the simulation study and plot the heatmap.
rm(list=ls())
library(scran)
library(cluster)

########################
# Load Simulation Data #
########################
# Working directory
# setwd("G:/scRNA/Journal/Github_reproduce/Simulation")

# Loading the file name list of all simulation count data
file_list <- list.files("./RawCountData")
file_list <- file_list[grepl("sim_count_data", file_list)]
B <- length(file_list)
SimulatedCounts <- list()
for(b in 1:B){
  file_name <- paste0("./RawCountData/",file_list[b])
  SimulatedCounts[[b]] <- as.matrix(read.table(file_name,header = T))
}

####################################
# Apply MNN to the Simulation Data #
####################################
data_MNN <- SimulatedCounts
# Referring to https://github.com/MarioniLab/MNN2017/blob/master/Pancreas
# Normalization
data_MNN_normalized<-list()
for(b in 1:B){
  high.abF <- scater::calcAverage(data_MNN[[b]]) > 1
  clustF <- quickCluster(data_MNN[[b]], min.size=10 , method="igraph", subset.row=high.abF)
  sizeF <- computeSumFactors(data_MNN[[b]], sizes=seq(11, 81, 5), cluster=clustF, subset.row=high.abF)
  data_MNN_normalized[[b]] <- t(t(data_MNN[[b]])/sizeF)
}

# Rescaling the first dataset to match the coverage of the second.
ave<-list()
for(b in 1:B){
  ave[[b]] <- rowMeans(data_MNN_normalized[[b]])
  if(b>1){
    data_MNN_normalized[[b]] <- data_MNN_normalized[[b]] * median(ave[[1]]/ave[[b]])
  }
}

# Performing log-transformation and save results to file.
log_data_MNN<-list()
for(b in 1:B){
  log_data_MNN[[b]] <- log1p(data_MNN_normalized[[b]])
}

# Performing the correction with MNN
mnn.out<-do.call(mnnCorrect,log_data_MNN)
X.mnn <- do.call(cbind,mnn.out$corrected)

# The corrected read count matrix
t.mnn <- t(X.mnn)

##############
# Clustering #
##############
# Referring to http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html
mnn.out.fast <- do.call(fastMNN, c(log_data_MNN, list(k=20, d=50, approximate=TRUE)))
dim(mnn.out.fast$corrected)
mnn.out.fast$batch@values <- paste0("Batch",1:B)
mnn.out.fast$pairs

time_spent <- end.time-start.time
print(time_spent)


omat <- do.call(cbind, log_data_MNN)
sce <- SingleCellExperiment(list(logcounts=omat))
reducedDim(sce, "MNN") <- mnn.out.fast$corrected
sce$Batch <- as.character(mnn.out.fast$batch)
sce


snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
clusters <- igraph::cluster_walktrap(snn.gr)
table(clusters$membership, sce$Batch)

# The estimated cell type indicators by MNN
w_MNN <- factor(clusters$membership)

# Store the workspace
if(!dir.exists("Workspace")){
  dir.create("Workspace")
}

save.image("./Workspace/MNN_workspace.RData")
