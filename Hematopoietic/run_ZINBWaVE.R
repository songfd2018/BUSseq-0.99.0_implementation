# Apply ZINB-WaVE to the hematopoietic study.
# Referring to https://github.com/drisso/zinb_analysis/blob/master/real_data/allen_covariates_1000.Rmd
rm(list=ls())
library(zinbwave)

set.seed(123)

###########################
# Load Hematopoietic Data #
###########################
# Working directory
# setwd("G:/scRNA/Journal/Github_reproduce/Mouse_Hematopoietic")

# Loading hematopoietic count data
load("./RawCountData/hemat_countdata.RData")
HematCounts <- list(GSE72857 = dataA2,
                    GSE81682 = dataF2)

B <- 2
nb <- c(ncol(dataA2),ncol(dataF2))

#############################################
# Apply ZINB-WaVE to the Hematopoietic Data #
#############################################

data_ZINBW <- NULL
for(b in 1:B){
  data_ZINBW <- cbind(data_ZINBW, HematCounts[[b]])
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