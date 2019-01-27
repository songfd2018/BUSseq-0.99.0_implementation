# Generate the heatmap of the estimated parameters and the logarithm of the estimated true read counts 
# as well as the logarithm of the corrected read counts
rm(list=ls())

library(BUSseq)
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

#########################
# Load BUSseq Workspace #
#########################
# Working directory
# setwd("G:/scRNA/Journal/Github_reproduce/Simulation")

load("./Workspace/BUSseq_workspace.RData")

######################################################################
# Draw the Heatmap for the True Parameter Values and True Count Data #
######################################################################
if(!dir.exists("Image")){
  dir.create("Image")
}

if(!dir.exists("Image/Heatmap")){
  dir.create("Image/Heatmap")
}

if(!dir.exists("Image/Other")){
  dir.create("Image/Other")
}


# Batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")

# Cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")

# Colors in the heatmap
colorsChoice <- colorRampPalette(c("#CFCFCF","black"))

##########################################################################
# The heatmap of the estimated log-scale mean expression levels (Fig 3f) #
log_mu_celltype_est <- log(celltype_mean_expression(BUSseqfits_simulation))

# Sorting cell-type-specific expression levels 
# in the decreasing order of the proprotions of the first cell type
# because the true cell types are sorted in the same decreasing order
cell_type_proportion_est <- BUSseqfits_simulation$pi.est
order.est <- order(cell_type_proportion_est[,1],decreasing = T)

# The splitting points for log-scale mean expression levels
break_logmusub <- seq(0, 6, 0.6)

png("./Image/Heatmap/heatmap_simulation_est_log_mean_expression_levels.png",width = 1080, height = 1440)
heatmap.3(log_mu_celltype_est[, order.est],
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_celltype,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(0.6,0.4,3.6),
          col=colorsChoice(10),breaks = break_logmusub,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"),keysize = 0.8, cexRow=0.5,trace = "none")
dev.off()

##################################################
# The heatmap of location batch effects (Fig 3g) #
batch_effects_est <- location_batch_effects(BUSseqfits_simulation)

# the splitting points for batch effects
break_nu <- seq(0, 6, 0.6)

png("./Image/Heatmap/heatmap_simulation_est_batch_effects.png",width = 1080, height = 1440)
heatmap.3(t(batch_effects_est),
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = color_by_batch,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key 
          lhei=c(0.6,0.4,3.6),
          col=colorsChoice(10),breaks = break_nu,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"),keysize = 0.8, cexRow=0.5,trace = "none")
dev.off()

#####################################################################
# The heatmap of log-scale underlying true read count data (Fig 3h) #
imputed_count_est <- imputed_read_counts(BUSseqfits_simulation)
log_imputed_count_est <- NULL
for(b in 1:B){
  log_imputed_count_est <- cbind(log_imputed_count_est, log1p(imputed_count_est[[b]]))
}

# Generating the double color bar
nb <- BUSseqfits_simulation$n.perbatch

# Batch color bar for count data matrix
color_by_batch_count <- rep(color_by_batch,nb)

# Cell type color bar for count data matrix
N <- BUSseqfits_simulation$n.cell
color_by_celltype_count <- rep(NA,N)
for(i in 1:N){
  color_by_celltype_count[i] <- color_by_celltype[which(w_BUSseq[i]==order.est)]
}


# Upper color bar for batch, lower color bar for cell type
col_annotation<-cbind(color_by_celltype_count,color_by_batch_count)
colnames(col_annotation) <- NULL

# The splitting points for log-scale count data
break_logcount <- seq(0,10,1)

png("./Image/Heatmap/heatmap_simulation_log_imputed_count.png",width = 1080, height = 1440)
heatmap.3(log_imputed_count_est,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.6,0.4,3.6),
          col=colorsChoice(10),breaks = break_logcount,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"),keysize = 0.8, cexRow=0.5,trace = "none")
dev.off()

###############################################################
# The heatmap of log-scale corrected read count data (Fig 3i) #

# Sharing the same double color bar and splitting points 
# as that of the underlying true count data

png("./Image/Heatmap/heatmap_simulation_log_corrected_count.png",width = 1080, height = 1440)
heatmap.3(log_corrected_count_est,
          dendrogram = "none",#with cluster tree
          Rowv = FALSE, Colv = FALSE,
          labRow = FALSE, labCol = FALSE,
          ColSideColors = col_annotation,
          lmat=rbind(c(5,4), c(0,1), c(3,2)),#1=heatmap, 2=row dendogram, 3=col dendogram, 4= key
          lhei=c(0.6,0.4,3.6),
          col=colorsChoice(10),breaks = break_logcount,
          density.info="histogram",
          hclustfun = function(c)hclust(c,method="average"),keysize = 0.8, cexRow=0.5,trace = "none")
dev.off()

####################################################################################################################
# The scatter plot of the estimated cell-specific size factors versus the true cell-specific size factors (Fig 3e) #

# Loading the true cell-specific size factors
cell_effect_true_vec <- unlist(read.table("./True_para/delta_syn.txt"))

# Obtaining the estimated cell-specific size factors
cell_effect_est <- cell_effect_values(BUSseqfits_simulation)
cell_effect_est_vec <- unlist(cell_effect_est)

png("./Image/Other/scatter_simulation_cell_size_factors.png",width = 540, height = 720)
par(mar = c(5.1,6.1,4.1,2.1)) 
plot(cell_effect_true_vec,cell_effect_est_vec,xlab = expression(paste("True ",delta)), ylab = expression(paste("Estimated ",delta)),type="n",ylim = c(-2,4),xlim = c(-2,4),cex.axis = 3, cex.lab = 3)
points(cell_effect_true_vec,cell_effect_est_vec,pch = 8,cex=3)
abline(a=0,b=1,lty=3,cex=3)
dev.off()