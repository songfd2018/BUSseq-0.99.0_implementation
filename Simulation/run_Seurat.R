# Apply Seurat to the simulation study.
# Referring to https://satijalab.org/seurat/immune_alignment.html
rm(list=ls())
library(Seurat)

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
G <- nrow(SimulatedCounts[[1]])


#######################################
# Apply Seurat to the Simulation Data #
#######################################
# Setting the Seurat Object
data_Seurat <- NULL
for(b in 1:B){
  data_Seurat <- cbind(data_Seurat, SimulatedCounts[[b]])
}
colnames(data_Seurat) <- paste0("Batch_",rep(1:B,nb),"_Cell_",c(1:nb[1],1:nb[2],1:nb[3],1:nb[4]),sep="")
rownames(data_Seurat) <- paste0("Gene_",1:G,sep="")
obj_Seurat <- CreateSeuratObject(raw.data = data_Seurat, project = "Simulation Study", min.cells = 5)

# QC and selecting cells for further analysis
obj_Seurat <- FilterCells(obj_Seurat, subset.names = "nGene", low.thresholds = 5, high.thresholds = Inf)

# Normalizing the data
obj_Seurat <- NormalizeData(obj_Seurat)

# Scaling the data
obj_Seurat <- ScaleData(obj_Seurat, display.progress = F)

# As multi-set Seurat requires more than two datasets, we will split our test object into
# four just for this example
Batch1 <- SubsetData(obj_Seurat,cells.use = obj_Seurat@cell.names[1:nb[1]])
Batch2 <- SubsetData(obj_Seurat,cells.use = obj_Seurat@cell.names[1:nb[2]+nb[1]])
Batch3 <- SubsetData(obj_Seurat,cells.use = obj_Seurat@cell.names[1:nb[3]+sum(nb[1:2])])
Batch4 <- SubsetData(obj_Seurat,cells.use = obj_Seurat@cell.names[1:nb[4]+sum(nb[1:3])])

Batch1@meta.data$group <- "Batch1"
Batch2@meta.data$group <- "Batch2"
Batch3@meta.data$group <- "Batch3"
Batch4@meta.data$group <- "Batch4"
obj_Seurat.list <- list(Batch1, Batch2, Batch3, Batch4)

# Performing a canonical correlation analysis on multiple batches
Seurat_result <- RunMultiCCA(object.list = obj_Seurat.list, genes.use = rownames(data_Seurat), num.ccs = 30)


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
                               dims.align = 1:15)

# Clustering
Seurat_result <- FindClusters(object = Seurat_result, reduction.type = "cca", dims.use = 1:15, 
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