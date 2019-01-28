# Apply Seurat to the Pancreas Study.
# Referring to https://satijalab.org/seurat/immune_alignment.html
rm(list=ls())
library(Seurat)

set.seed(12345)

######################
# Load Pancreas Data #
######################
# Working directory
# setwd("G:/scRNA/Journal/Github_reproduce/Human_Pancreas")


# Loading the file name list of all simulation count data
load("./RawCountData/pancreas_countdata.RData")


B <- length(PancreasCounts)
nb <- rep(NA,B)
for(b in 1:B){
  nb[b] <- ncol(PancreasCounts[[b]])
}
G <- nrow(PancreasCounts[[1]])

##########################################
# Apply Seurat to the Hematopoietic Data #
##########################################
# Setting the Seurat Object
data_Seurat <- NULL
for(b in 1:B){
  data_Seurat <- cbind(data_Seurat, PancreasCounts[[b]])
}

obj_Seurat <- CreateSeuratObject(raw.data = data_Seurat, project = "Pancreas Study", min.cells = 5)

# QC and selecting cells for further analysis
obj_Seurat <- FilterCells(obj_Seurat, subset.names = "nGene", low.thresholds = 5, high.thresholds = Inf)

# Normalizing the data
obj_Seurat <- NormalizeData(obj_Seurat)

# Scaling the data
obj_Seurat <- ScaleData(obj_Seurat, display.progress = F)

# As multi-set Seurat requires more than two datasets, we will split our test object into
# four just for this example
Batch_data <- list()
for(b in 1:B){
  if(b==1){
    Batch_data[[b]] <- SubsetData(obj_Seurat,cells.use = obj_Seurat@cell.names[1:nb[1]])
  }else{
    Batch_data[[b]] <- SubsetData(obj_Seurat,cells.use = obj_Seurat@cell.names[1:nb[b] + sum(nb[1:(b-1)])])
  }
  print(Batch_data[[b]])
  Batch_data[[b]]@meta.data$group <- paste("Batch",b,sep = "_")
}

# Performing a canonical correlation analysis on multiple batches
Seurat_result <- RunMultiCCA(object.list = Batch_data, genes.use = rownames(data_Seurat), num.ccs = 30)

# Visualizing results of Seurat plot CC1 versus CC2 and look at a violin plot
p1 <- DimPlot(object = Seurat_result, reduction.use = "cca", group.by = "group", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = Seurat_result, features.plot = "CC1", group.by = "group", 
              do.return = TRUE)
plot_grid(p1, p2)

# Choosing CCs for downstream analysis
p3 <- MetageneBicorPlot(Seurat_result, grouping.var = "group", dims.eval = 1:30, 
                        display.progress = FALSE)


# Aligning the CCA subspaces
Seurat_result <- AlignSubspace(Seurat_result, reduction.type = "cca", grouping.var = "group", 
                               dims.align = 1:12)

# Clustering
Seurat_result <- FindClusters(object = Seurat_result, reduction.type = "cca", dims.use = 1:12, 
                              resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = Seurat_result)

w_Seurat<-Seurat_result@meta.data$res.0.6

# Extracting the first 15 canonical correlations
extract.cca <- GetCellEmbeddings(object = Seurat_result, reduction.type = "cca.aligned")

# Store the workspace
if(!dir.exists("Workspace")){
  dir.create("Workspace")
}

save.image("./Workspace/Seurat_workspace.RData")