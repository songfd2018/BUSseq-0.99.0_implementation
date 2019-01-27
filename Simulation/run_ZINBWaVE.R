# Apply ZINB-WaVE to the simulation study and plot the heatmap.
# Referring to https://github.com/drisso/zinb_analysis/blob/master/real_data/allen_covariates_1000.Rmd
rm(list=ls())
library(zinbwave)

set.seed(12345)

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
nb <- rep(NA,B)
for(b in 1:B){
  file_name <- paste0("./RawCountData/",file_list[b])
  SimulatedCounts[[b]] <- as.matrix(read.table(file_name,header = T))
  nb[b] <- ncol(SimulatedCounts[[b]])
}



##########################################
# Apply ZINB-WaVE to the Simulation Data #
##########################################

data_ZINBW <- NULL
for(b in 1:B){
  data_ZINBW <- cbind(data_ZINBW, SimulatedCounts[[b]])
}

# Factorizing the batch indicators
batch_ind <- factor(rep(1:B,nb))

# Performing ZINB-WaVE
zinb_batch <- zinbFit(data_ZINBW, K = 10, X=model.matrix(~batch_ind), epsilon=1e3)

# Clustering
library("clusterExperiment")

ZINBW_clust<-clusterSingle(t(zinb_batch@W),sequential= T, subsample =F,mainClusterArgs=list(clusterFunction="kmeans"),seqArgs=list(k0=5,beta=0.95))
w_ZINBW<-ZINBW_clust@clusterMatrix

# Store the workspace
if(!dir.exists("Workspace")){
  dir.create("Workspace")
}

save.image("./Workspace/ZINBWaVE_workspace.RData")