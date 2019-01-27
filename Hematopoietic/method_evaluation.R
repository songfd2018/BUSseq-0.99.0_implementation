rm(list=ls())
library(mclust)
library(ggplot2)
library(Rtsne)

###################
# Load Workspaces #
###################
# Working directory
# setwd("G:/scRNA/Journal/Github_reproduce/Mouse_Hematopoietic")

load("./Workspace/BUSseq_workspace.RData")
load("./Workspace/MNN_workspace.RData")
load("./Workspace/Seurat_workspace.RData")
load("./Workspace/ZINBWaVE_workspace.RData")

#############################################################
# Evaluation of the Performance of All the Methods          #
# according to ARI, silhouette coefficients and t-SNE plots #
#############################################################
if(!dir.exists("Results")){
  dir.create("Results")
}

#######
# ARI #
#######
ARI_values <- rep(NA, 4)
names(ARI_values) <- c("BUSseq", "MNN", "Seurat", "ZINB-WaVE")

# Loading the true cell type indicators
celltype_true <- NULL
for(b in 1:B){
  celltype_true <- c(celltype_true, colnames(HematCounts[[b]]))
}
N <- length(celltype_true)
w_true <- rep(NA, N)
w_true[which(celltype_true=="CMP")]<-1
w_true[which(celltype_true=="ERY" | celltype_true == "MEP")]<-2
w_true[which(celltype_true=="GMP")]<-3
w_true[which(celltype_true=="LMPP")]<-4
w_true[which(celltype_true=="LTHSC")]<-5
w_true[which(celltype_true=="MPP")]<-6
w_true[which(celltype_true=="other")]<-7

# Calculating the ARI
ARI_values[1] <- adjustedRandIndex(w_true,w_BUSseq)
ARI_values[2] <- adjustedRandIndex(w_true,w_MNN)
ARI_values[3] <- adjustedRandIndex(w_true,w_Seurat)
ARI_values[4] <- adjustedRandIndex(w_true,w_ZINBW)

write.table(ARI_values, "./Results/ARI_values.txt",col.names = F)

###########################
# Silhouette Coefficients #
###########################
silhouette_cal<-function(Reads,#N * G
                         sub_ind){
  
  N<- nrow(Reads)
  sil_coef <- rep(NA,N)
  all.dists <- as.matrix(dist(Reads))
  
  for(i in 1:N){
    clust_index <- which(sub_ind==sub_ind[i])
    other_index <- which(sub_ind!=sub_ind[i])
    aver_dist <- sum(all.dists[i,clust_index])/length(clust_index-1)
    min_dist <- min(all.dists[i,other_index])
    sil_coef[i] <-  (min_dist - aver_dist)/max(aver_dist,min_dist)
  }
  return(sil_coef)
}

# Silhouette coefficients of the corrected read count matrix of 
# the identified intrinsic genes by BUSseq
sil_BUSseq <- silhouette_cal(t(log_corrected_count_est[intrinsic_gene_indices,]),w_BUSseq)

# Silhouette coefficients of the whole corrected read count matrix from MNN
sil_MNN <- silhouette_cal(t.mnn, w_MNN)

# Silhouette coefficients of the top 25 CCs inferred by Seurat
sil_Seurat <- silhouette_cal(extract.cca, w_Seurat)

# Silhouette coefficients of the top 10-component W inferred by ZINB-WaVE
sil_ZINBW <- silhouette_cal(zinb_batch@W, w_ZINBW)

# Storing the Silhouette coefficients
sil_matrix <- cbind(sil_BUSseq,sil_MNN,sil_Seurat,sil_ZINBW)
colnames(sil_matrix) <- c("BUSseq","MNN","Seurat","ZINB-WaVE")
write.table(sil_matrix,"./Results/Silhouette_coef.txt",row.names = F)

# Drawing the violin plot
sil_com<-data.frame(sli_coef=c(sil_BUSseq,sil_MNN,sil_ZINBW,sil_Seurat),
                    method=factor(rep(c("BUSseq","MNN","ZINB-WaVE","Seurat"),each=N)))

png("./Image/Other/violin_plot_of_Silhouette_coef.png",width = 960, height = 720)
p <- ggplot(sil_com, aes(x=method, y=sli_coef,fill=method)) + 
  geom_violin() + xlab(NULL) + ylab(NULL)+theme_bw() +
  theme(text = element_text(size=32),axis.text.x =element_text(size=24), axis.text.y = element_text(size =24),panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
p
dev.off()

##############
# t-SNE plot #
##############
if(!dir.exists("Image/tSNE")){
  dir.create("Image/tSNE")
}

plot_tSNE<-function(tSNE_var,col,p_name,leg_name){
  
  # The image with legend
  png(paste(p_name,"_with_legend.png",sep=""),width = 1440, height = 1080)
  par(mar=c(5.1,6.1,4.1,2.1))
  plot(tSNE_var,t="n",xaxt="n",yaxt="n",xlab="",ylab="")#main="tsne_uncorrected_by_batch")
  axis(1,cex.axis=6,line=2.5,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  points(tSNE_var,pch=19 ,col= col,cex=3)
  legend("bottomright",legend=leg_name[,1],pch=19,col=leg_name[,2],cex=2)
  dev.off()
  
  # The image without legend
  png(paste(p_name,".png",sep=""),width = 1440, height = 1080)
  par(mar=c(5.1,6.1,4.1,2.1))
  plot(tSNE_var,t="n",xaxt="n",yaxt="n",xlab="",ylab="")#main="tsne_uncorrected_by_batch")
  axis(1,cex.axis=6,line=2.5,tick = F)#plot the x axis
  axis(2,cex.axis=6,tick = F)#plot the y axis
  points(tSNE_var,pch=19 ,col= col,cex=3)
  dev.off()
  
}

#############
#Uncorrected#
#############
data_raw <- NULL
for(b in 1:B){
  data_raw <- cbind(data_raw, HematCounts[[b]])
}

set.seed(123)
all.dists.unc <- as.matrix(dist(log1p(t(data_raw))))
tsne_uncorrected <- Rtsne(all.dists.unc, is_distance=TRUE, perplexity = 30)

########
#BUSseq#
########
set.seed(123)
all.dists.BUSseq <- as.matrix(dist(t(log_corrected_count_est[intrinsic_gene_indices,])))
tsne_BUSseq_dist <- Rtsne(all.dists.BUSseq, is_distance=TRUE, perplexity = 30)

#####
#MNN#
#####
set.seed(123)
all.dists.mnn <- as.matrix(dist(t.mnn))
tsne_MNN_dist <- Rtsne(all.dists.mnn, is_distance=TRUE, perplexity = 30)

########
#Seurat#
########
set.seed(123)
all.dists.Seurat <- as.matrix(dist(extract.cca))
tsne_Seurat <- Rtsne(all.dists.Seurat, is_distance=TRUE, perplexity = 30)

###########
#ZINB-WaVE#
###########
set.seed(123)
all.dists.ZINBW <- as.matrix(dist(zinb_batch@W))
tsne_ZINBW <- Rtsne(all.dists.ZINBW, is_distance=TRUE, perplexity = 30)

###################################################
# T-SNE plot colored by batch and true cell types #
###################################################
# Batch color bar
color_by_batch<-c("#EB4334","#4586F3")
# Transparenting the color
color_by_batch_trans<-paste0(color_by_batch,"44",sep="")

# Cell type color bar
# Interpolating the set of cell type colors to represent more cell types 
color_by_celltype_ramp<-colorRampPalette(c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE"))(19)
# Notice that color_by_celltype is the same as color_by_celltype_ramp[c(1,4,7,10,13,16,19)]
color_by_celltype<-color_by_celltype_ramp[c(1,4,7,10,13,16,19)]


cell_color_by_batch <-  color_by_batch_trans[rep(1:B,nb)]
cell_color_by_celltype <- color_by_celltype[w_true]

# legend information 
batch_name <- c("GSE72857","GSE81682")
celltype_name <- c("CMP","MEP","GMP","LMPP","LTHSC","MPP","other")
legend_batch <- cbind(batch_name,color_by_batch)
legend_celltype <- cbind(celltype_name,color_by_celltype)

#############
#Uncorrected#
#############
unc_by_batch <- "./Image/tSNE/tsne_simulation_uncorrected_by_batch"
plot_tSNE(tsne_uncorrected$Y,cell_color_by_batch,unc_by_batch,legend_batch)

unc_by_celltype<- "./Image/tSNE/tsne_simulation_uncorrected_by_celltype"
plot_tSNE(tsne_uncorrected$Y,cell_color_by_celltype,unc_by_celltype,legend_celltype)

########
#BUSseq#
########
BUSseq_by_batch <- "./Image/tSNE/tsne_simulation_BUSseq_by_batch"
plot_tSNE(tsne_BUSseq_dist$Y,cell_color_by_batch,BUSseq_by_batch,legend_batch)

BUSseq_by_celltype<- "./Image/tSNE/tsne_simulation_BUSseq_by_celltype"
plot_tSNE(tsne_BUSseq_dist$Y,cell_color_by_celltype,BUSseq_by_celltype,legend_celltype)

#####
#MNN#
#####
MNN_by_batch <- "./Image/tSNE/tsne_simulation_MNN_by_batch"
plot_tSNE(tsne_MNN_dist$Y,cell_color_by_batch,MNN_by_batch,legend_batch)

MNN_by_celltype<- "./Image/tSNE/tsne_simulation_MNN_by_celltype"
plot_tSNE(tsne_MNN_dist$Y,cell_color_by_celltype,MNN_by_celltype,legend_celltype)

########
#Seurat#
########
Seurat_by_batch <-  "./Image/tSNE/tsne_simulation_Seurat_by_batch"
plot_tSNE(tsne_Seurat$Y,cell_color_by_batch,Seurat_by_batch,legend_batch)

Seurat_by_celltype<- "./Image/tSNE/tsne_simulation_Seurat_by_celltype"
plot_tSNE(tsne_Seurat$Y,cell_color_by_celltype,Seurat_by_celltype,legend_celltype)

###########
#ZINB-WaVE#
###########
ZINBW_by_batch <- "./Image/tSNE/tsne_simulation_ZINBWaVE_by_batch"
plot_tSNE(tsne_ZINBW$Y,cell_color_by_batch,ZINBW_by_batch,legend_batch)

ZINBW_by_celltype<- "./Image/tSNE/tsne_simulation_ZINBWaVE_by_celltype"
plot_tSNE(tsne_ZINBW$Y,cell_color_by_celltype,ZINBW_by_celltype,legend_celltype)

##############################################
# T-SNE plot colored by estimated cell types #
##############################################
# Adjusting the order of the estimated cell types to obtain a better correspondence
# between the true cell type indicators and the estimated cell type indicators
table(w_true, w_BUSseq)
w_BUSseq_shift <- w_BUSseq
w_BUSseq_shift[which(w_BUSseq=="1")]<-1#1 in w
w_BUSseq_shift[which(w_BUSseq=="4")]<-2#2 in w
w_BUSseq_shift[which(w_BUSseq=="2")]<-3#3 in w
w_BUSseq_shift[which(w_BUSseq=="3")]<-4#4 in w
w_BUSseq_shift[which(w_BUSseq=="7")]<-5#5 in w
w_BUSseq_shift[which(w_BUSseq=="5")]<-6#6 in w
w_BUSseq_shift[which(w_BUSseq=="6")]<-7#7 in w
BUSseq_col <- color_by_celltype_ramp[c(1,4,7,10,13,16,19)]

table(w_true,w_MNN)
w_MNN_shift<-w_MNN
w_MNN_shift[which(w_MNN=="3")]<-1#1 in w
w_MNN_shift[which(w_MNN=="5")]<-2
w_MNN_shift[which(w_MNN=="9")]<-3
w_MNN_shift[which(w_MNN=="2")]<-4#2 in w
w_MNN_shift[which(w_MNN=="8")]<-5
w_MNN_shift[which(w_MNN=="4")]<-6
w_MNN_shift[which(w_MNN=="7")]<-7#3 in w
w_MNN_shift[which(w_MNN=="6")]<-8
w_MNN_shift[which(w_MNN=="13")]<-9
w_MNN_shift[which(w_MNN=="15")]<-10#4 in w
w_MNN_shift[which(w_MNN=="14")]<-11
w_MNN_shift[which(w_MNN=="11")]<-12
w_MNN_shift[which(w_MNN=="12")]<-13#5 in w
w_MNN_shift[which(w_MNN=="16")]<-14#6 in w
w_MNN_shift[which(w_MNN=="10")]<-15
w_MNN_shift[which(w_MNN=="1")]<-16#7 in w
MNN_col <- color_by_celltype_ramp[c(1,2,3,4,5,6,7,8,9,10,11,12,13,16,18,19)]

table(w_true,w_Seurat)
w_Seurat_shift<-w_Seurat
w_Seurat_shift[which(w_Seurat=="3")]<-1#1 in w
w_Seurat_shift[which(w_Seurat=="5")]<-2
w_Seurat_shift[which(w_Seurat=="9")]<-3
w_Seurat_shift[which(w_Seurat=="0")]<-4#2 in w
w_Seurat_shift[which(w_Seurat=="7")]<-5
w_Seurat_shift[which(w_Seurat=="2")]<-6#3 in w
w_Seurat_shift[which(w_Seurat=="4")]<-7
w_Seurat_shift[which(w_Seurat=="8")]<-8#4 in w
w_Seurat_shift[which(w_Seurat=="10")]<-9#5 in w
w_Seurat_shift[which(w_Seurat=="1")]<-10#6 in w
w_Seurat_shift[which(w_Seurat=="6")]<-11#7 in w
Seurat_col <- color_by_celltype_ramp[c(1,2,3,4,6,7,8,10,13,16,19)]

table(w_true,w_ZINBW)
w_ZINBW_shift<-w_ZINBW
w_ZINBW_shift[which(w_ZINBW=="-1")]<-1#1 in w
w_ZINBW_shift[which(w_ZINBW=="10")]<-2
w_ZINBW_shift[which(w_ZINBW=="11")]<-3
w_ZINBW_shift[which(w_ZINBW=="1")]<-4#2 in w_true
w_ZINBW_shift[which(w_ZINBW=="3")]<-5
w_ZINBW_shift[which(w_ZINBW=="5")]<-6
w_ZINBW_shift[which(w_ZINBW=="2")]<-7#3 in w
w_ZINBW_shift[which(w_ZINBW=="4")]<-8#4 in w_true
w_ZINBW_shift[which(w_ZINBW=="7")]<-9#5 in w_true
w_ZINBW_shift[which(w_ZINBW=="9")]<-10
w_ZINBW_shift[which(w_ZINBW=="8")]<-11#6 in w_true
w_ZINBW_shift[which(w_ZINBW=="6")]<-12#7 in w_true
ZINBW_col <- color_by_celltype_ramp[c(1,2,3,4,5,6,7,10,13,15,16,19)]

#BUSseq
col_by_celltype_est <- BUSseq_col[w_BUSseq_shift]
BUSseq_by_celltype_est<- "./Image/tSNE/tsne_simulation_BUSseq_by_celltype_est"
legend_celltype_BUSseq<- cbind(paste0("Cell type ",1:length(table(w_BUSseq_shift)),sep=""),BUSseq_col)
plot_tSNE(tsne_BUSseq_dist$Y,col_by_celltype_est,BUSseq_by_celltype_est,legend_celltype_BUSseq)

#MNN
col_by_celltype_est <-  MNN_col[w_MNN_shift]
MNN_by_celltype_est<- "./Image/tSNE/tsne_simulation_MNN_by_celltype_est"
legend_celltype_MNN<- cbind(paste0("Cell type ",1:length(table(w_MNN_shift)),sep=""),MNN_col)
plot_tSNE(tsne_MNN_dist$Y,col_by_celltype_est,MNN_by_celltype_est,legend_celltype_MNN)

#Seurat
col_by_celltype_est <- Seurat_col[as.numeric(w_Seurat_shift)]
Seurat_by_celltype_est<- "./Image/tSNE/tsne_simulation_Seurat_by_celltype_est"
legend_celltype_Seurat<- cbind(paste0("Cell type ",1:length(table(w_Seurat_shift)),sep=""),Seurat_col)
plot_tSNE(tsne_Seurat$Y,col_by_celltype_est,Seurat_by_celltype_est,legend_celltype_Seurat)

#ZINB-WaVE
col_by_celltype_est <- ZINBW_col[as.numeric(w_ZINBW_shift)]
ZINBW_by_celltype_est<- "./Image/tSNE/tsne_simulation_ZINBWaVE_by_celltype_est"
legend_celltype_ZINBW<- cbind(paste0("Cell type ",1:length(table(w_ZINBW_shift)),sep=""),ZINBW_col)
plot_tSNE(tsne_ZINBW$Y,col_by_celltype_est,ZINBW_by_celltype_est,legend_celltype_ZINBW)

# Store the workspace
if(!dir.exists("Workspace")){
  dir.create("Workspace")
}

save.image("./Workspace/comparison_workspace.RData")