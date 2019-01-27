rm(list=ls())
library("devtools")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# setwd("G:/scRNA/Journal/Github_reproduce/Simulation") # temp

set.seed(12356)
################################
# Set the Synthetic Parameters #
################################
#The number of batches
B<-4

#The number of cells per batch
nb<-c(150,150,150,150)

#The total number of cells
N<-sum(nb)

#The number of genes
G<-10000

#The number of cell types
K<-5

#The first column of gamma.syn denotes the intercept of 
#the logistic regression for dropout events
#The second column of gamma.syn denotes the odds ratios 
#of the logistic regression for dropout events
gamma.syn<-matrix(0,B,2)
gamma.syn[1,]<-c(1.5,0.01)
gamma.syn[2,]<-c(1.5,0.01)
gamma.syn[3,]<-c(1,0.02)
gamma.syn[4,]<-c(1,0.02)

#the log-scale baseline expression levels
alpha.syn<-rep(1,G)
alpha.syn[1:150]<-3

#the cell-type effects 
beta.syn<-matrix(NA,G,K)

#the first cell type is regarded as the reference cell type 
#without cell-type effects 
beta.syn[,1] <- 0

#the cell-type effects of the second cell type
beta.syn[1:100,2] <- -2
beta.syn[101:150,2] <- 2
beta.syn[151:250,2] <- 3
beta.syn[251:G,2] <- 0

#the cell-type effects of the third cell type
beta.syn[1:150,3]<- -2
beta.syn[151:170,3] <- 0
beta.syn[171:300,3] <- 2
beta.syn[301:G,3] <- 0

#the cell-type effects of the forth cell type
beta.syn[1:150,4]<- -2
beta.syn[151:250,4] <- 0
beta.syn[251:320,4] <- 3
beta.syn[321:400,4] <- 2.5
beta.syn[401:G,4] <- 0

#the cell-type effects of the fifth cell type
beta.syn[1:350,5] <- 0
beta.syn[351:500,5] <- 4
beta.syn[501:G,5] <- 0

#the batch effects
nu.syn<-matrix(NA,B,G)

#the first batch is taken as the reference batch 
#without batch effects
nu.syn[1,] <- 0

#the batch effect of the second batch
nu.syn[2,1:4000] <- 3
nu.syn[2,4001:8000] <- 1
nu.syn[2,6001:10000] <- 3

#the batch effect of the third batch
nu.syn[3,1:3000] <- 3
nu.syn[3,3001:6000] <- 2
nu.syn[3,6001:10000] <- 1

#the batch effect of the forth batch
nu.syn[4,1:4000] <- 2
nu.syn[4,4001:8000] <-3
nu.syn[4,8001:10000] <- 4

#the cell-specific size factors
delta.syn <- list()
for(b in 1:B){
  delta.syn[[b]] <- rep(NA, nb[b])
}


#the first cell in each batch is regarded as the reference cell 
#with the cell-specific size factors being 0
delta.syn[[1]][1:50] <- 0
delta.syn[[1]][51:90] <- 1
delta.syn[[1]][91:150] <- 2

#the second batch
delta.syn[[2]][1:50] <- 0
delta.syn[[2]][51:70] <- 2
delta.syn[[2]][71:100] <- -1
delta.syn[[2]][101:150] <- 1

#the third batch
delta.syn[[3]][1:30] <- 0
delta.syn[[3]][31:90] <- 1
delta.syn[[3]][91:150] <- 2

#the forth batch
delta.syn[[4]][1:70] <- 0
delta.syn[[4]][71:150] <- 3

#the batch-specific and gene-specific overdispersion parameters
phi.syn<-matrix(10,B,G)#mean 2 var 0.5
phi.syn[,8001:10000]<-1

#the cell-type proportions in each batch
pi.syn <- matrix(NA,K,B)

#the first batch
pi.syn[,1]<-c(0.4,0.3,0.2,0.1,0)

#the second batch
pi.syn[,2]<-c(0,0.2,0.3,0.3,0.2)

#the third batch
pi.syn[,3]<-c(0.24,0,0.2,0.26,0.3)

#the forth batch
pi.syn[,4]<-c(0,0.3,0.4,0.3,0)


##############################################
# Simulate Latent Varibles and Observed data #
##############################################
#the cell-type indicators of each cell
w <- list()

#the first batch
w[[1]] <- rep(1:K, nb[1] * pi.syn[,1])

#the second batch
w[[2]] <- rep(1:K, nb[2] * pi.syn[,2])

#the third batch
w[[3]] <- rep(1:K, nb[3] * pi.syn[,3])

#the forth batch
w[[4]] <- rep(1:K, nb[4] * pi.syn[,4])


#the indicators for dropout events
z<-list()

#the underlying true expression levels
x<-list()

#the observed expression levels
y<-list()

#the logarithm of mean expreesion level of each gene in each cell
log.mu<-list()

for(b in 1:B){
  z[[b]] <- matrix(NA, G, nb[b])
  x[[b]] <- matrix(NA, G, nb[b])
  y[[b]] <- matrix(NA, G, nb[b])
  log.mu[[b]] <- matrix(NA, G, nb[b])
}

#generate the latent variable and observed data
for(b in 1:B){
  for(i in 1:nb[b]){
    log.mu[[b]][,i] <- alpha.syn + beta.syn[,w[[b]][i]] 
    log.mu[[b]][,i] <- log.mu[[b]][,i] + nu.syn[b,]
    log.mu[[b]][,i] <- log.mu[[b]][,i] + delta.syn[[b]][i]
    
    for(j in 1:G){
      x[[b]][j,i]<-rnbinom(1,phi.syn[b,j],
                           mu=exp(log.mu[[b]][j,i]))
      logit_pi <- gamma.syn[b,1] + gamma.syn[b,2] * x[[b]][j,i]
      
      z[[b]][j,i]<-rbinom(1,1,prob = 1/(1+exp(-logit_pi)))
      if(z[[b]][j,i]==1){
        y[[b]][j,i]<- x[[b]][j,i]
      }else{
        y[[b]][j,i]<- 0
      }
    }
  }
}

######################################################################
# Draw the Heatmap for the True Parameter Values and True Count Data #
######################################################################
if(!dir.exists("Image")){
  dir.create("Image")
}

if(!dir.exists("Image/Heatmap")){
  dir.create("Image/Heatmap")
}

# batch color bar
color_by_batch<-c("#EB4334","#FBBD06","#35AA53","#4586F3")

# cell type color bar
color_by_celltype<-c("#F88A7E", "#FFD87D", "#ABD978","#8097D3","#9C7ACE")

# colors in the heatmap
colorsChoice <- colorRampPalette(c("#CFCFCF","black"))

############################################################
# the heatmap of log-scale mean expression levels (Fig 3a) #
log_mu_celltype <- alpha.syn + beta.syn

# the splitting points for log-scale mean expression levels
break_logmusub <- seq(0, 6, 0.6)

png("./Image/Heatmap/heatmap_simulation_log_mean_expression_levels.png",width = 1080, height = 1440)
heatmap.3(log_mu_celltype,
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
# the heatmap of location batch effects (Fig 3b) #
# the splitting points for batch effects
break_nu <- seq(0, 6, 0.6)

png("./Image/Heatmap/heatmap_simulation_batch_effects.png",width = 1080, height = 1440)
heatmap.3(t(nu.syn),
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
# the heatmap of log-scale underlying true read count data (Fig 3c) #
log_x <- NULL
for(b in 1:B){
  log_x <- cbind(log_x, log1p(x[[b]]))
}

#generate the double color bar
#batch color bar for count data matrix
color_by_batch_count <- rep(color_by_batch,nb)

#cell type color bar for count data matrix
color_by_celltype_count <- NULL
for(b in 1:B){
  color_by_celltype_count <- c(color_by_celltype_count, color_by_celltype[w[[b]]])
}

#upper color bar for batch, lower color bar for cell type
col_annotation<-cbind(color_by_celltype_count,color_by_batch_count)
colnames(col_annotation) <- NULL

#the splitting points for log-scale count data
break_logcount <- seq(0,10,1)

png("./Image/Heatmap/heatmap_simulation_log_underlying_true_count.png",width = 1080, height = 1440)
heatmap.3(log_x,
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

##############################################################
# the heatmap of log-scale observed read count data (Fig 3d) #
log_y <- NULL
for(b in 1:B){
  log_y <- cbind(log_y, log1p(y[[b]]))
}

#share the same double color bar and splitting points 
#as that of the underlying true count data

png("./Image/Heatmap/heatmap_simulation_log_observed_count.png",width = 1080, height = 1440)
heatmap.3(log_y,
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

##################################
# Save the Simulation Count Data #
##################################
if(!dir.exists("RawCountData")){
  dir.create("RawCountData")
}

for(b in 1:B){
  file_name <- paste0("./RawCountData/sim_count_data_batch",b,".txt")
  write.table(y[[b]],file_name)
}

##################################################
# Save the Cell-specific Size Factors for Fig 3e #
##################################################
if(!dir.exists("True_para")){
  dir.create("True_para")
}

#write delta.syn as a column vector
delta.file <- "./True_para/delta_syn.txt"
file.create(delta.file)

for(b in 1:B){
  write.table(delta.syn[[b]],delta.file,append = T, col.names = F, row.names = F)
}

#write w as a column vector
w.file <- "./True_para/w_syn.txt"
file.create(w.file)


for(b in 1:B){
  write.table(w[[b]],w.file,append = T, col.names = F, row.names = F)
}