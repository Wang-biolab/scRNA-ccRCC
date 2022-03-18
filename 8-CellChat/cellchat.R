library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(PCAtools)
options(stringsAsFactors = FALSE)
rm(list=ls())

load('scRNA_celltype.Rdata')
table(scRNA$orig, scRNA$celltype)
##remove platelet
Idents(scRNA) <- 'celltype'
scRNA <- subset(scRNA, idents = c(
  'CD14','CD4+T','CD8+T_1','B','NK','CD16','cDC2','CD8+T_2',
  'NKT','Treg','CD8+_T3', 'pDC', 'cDC1',
  'Plasmablast', 'HSPC'
))
scRNA$celltype <- as.factor(as.character(scRNA$celltype))
table(scRNA$orig)
Idents(scRNA) <- 'orig'

tumor <- subset(scRNA, idents = 'T')
normal <- subset(scRNA, idents = 'N')
pbmc <- subset(scRNA, idents = 'P')


### 创建cellchat对象
cco.tumor <- createCellChat(tumor@assays$RNA@data, meta = tumor@meta.data, group.by = 'celltype')
cco.normal <- createCellChat(normal@assays$RNA@data, meta = normal@meta.data, group.by = 'celltype')
cco.pbmc <- createCellChat(pbmc@assays$RNA@data, meta = pbmc@meta.data, group.by = 'celltype')
save(cco.tumor, cco.normal, cco.pbmc, file ='cco_T_N_P.RData')

## Part1 ###
cellchat <- cco.tumor
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cell = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')
cco.tumor<- cellchat
save(cco.tumor, file = 'cco.tumor.after_process.Rdata')

#netAnalysis_signalingRole_network(cco.pbmc, signaling = 'MIF', width = 13, height = 4, font.size = 13)
## Part2 ###
cellchat <- cco.normal
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cell = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')
cco.normal <- cellchat
save(cco.normal, file = 'cco.normal.after_process.Rdata')

cellchat <- cco.pbmc
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
#cellchat <- filterCommunication(cellchat, min.cell = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')
cco.pbmc<- cellchat
save(cco.pbmc, file = 'cco.pbmc.after_process.Rdata')

cco.list <- list( Normal = cco.normal,Tumor = cco.tumor,PBMC = cco.pbmc)

cellchat <- mergeCellChat(cco.list, add.names = names(cco.list),
                          cell.prefix = TRUE)


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),
                           measure = 'count'
)

gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2),
                           measure = 'weight'
)
p <- gg1 + gg2
p
ggsave('Overview_number_strength.png',
       p, width = 6, height = 3.5)


par(mfrow = c(1,1))
groupSize <- as.numeric(table(cco.tumor@idents))
netVisual_circle(cco.tumor@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title = 'Number of interactions in tumor'
)
groupSize <- as.numeric(table(cco.normal@idents))
netVisual_circle(cco.normal@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title = 'Number of interactions in normal'
)

groupSize <- as.numeric(table(cco.pbmc@idents))
netVisual_circle(cco.pbmc@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title = 'Number of interactions in PBMC'
)

groupSize <- as.numeric(table(cco.tumor@idents))
netVisual_circle(cco.tumor@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title = 'Strength of interactions in tumor'
)
groupSize <- as.numeric(table(cco.normal@idents))
netVisual_circle(cco.normal@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title = 'Strength of interactions in normal'
)

groupSize <- as.numeric(table(cco.pbmc@idents))
netVisual_circle(cco.pbmc@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F,
                 title = 'Strength of interactions in PBMC'
)


############ Net strength individual in normal #####
mat <- cco.normal@net$weight
groupSize <- as.numeric(table(cco.normal@idents))
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat),
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   arrow.width = 0.2,
                   arrow.size = 0.1,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
par("mar")
par(mar = c(1,1,1,1))
######################################

############ Net strength individual in tumor #####
mat <- cco.tumor@net$weight
groupSize <- as.numeric(table(cco.tumor@idents))
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat),
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   arrow.width = 0.2,
                   arrow.size = 0.1,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
par("mar")
par(mar = c(1,1,1,1))
######################################
############ Net strength individual in PBMC #####
mat <- cco.pbmc@net$weight
groupSize <- as.numeric(table(cco.pbmc@idents))
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat),
                 ncol = ncol(mat),
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = T,
                   arrow.width = 0.2,
                   arrow.size = 0.1,
                   edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}
par("mar")
par(mar = c(1,1,1,1))
################# Merge Normal vs. Tumor ##################

cco.list <- list( Normal = cco.normal,Tumor = cco.tumor)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list),
                          cell.prefix = TRUE)
par(mfrow = c(1,1), xpd = TRUE)
netVisual_diffInteraction(cellchat,
                          weight.scale = T,
                          measure = 'weight',
                          title.name = 'Differential interaction strength between N v.s.T'
)

h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = 'weight',
                        title.name = 'interaction heatmap between N v.s.T')
h1 + h2

num.link <- sapply(cco.list, function(x){
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) #control the dot size in the different datasets
gg <- list()
for(i in 1:length(cco.list)){
  gg[[i]] <- netAnalysis_signalingRole_scatter(cco.list[[i]],
                                               title = names(cco.list)[i],
                                               weight.MinMax = weight.MinMax)}
patchwork::wrap_plots(plots = gg)


cellchat <- computeNetSimilarityPairwise(cellchat, type = 'functional')
cellchat <- netEmbedding(cellchat, type = 'functional')
cellchat <- netClustering(cellchat, type = 'functional')
netVisual_embeddingPairwise(cellchat,
                            type = 'functional',
                            label.size = 3.5,
                            title = 'Compute signaling network similarity for N and T'
)


cellchat <- computeNetSimilarityPairwise(cellchat, type = 'structural')
cellchat <- netEmbedding(cellchat, type = 'structural')
cellchat <- netClustering(cellchat, type = 'structural')
netVisual_embeddingPairwise(cellchat, type = 'structural', label.size = 3.5,
                            title = 'Compute structural similarity for N and T')
rankSimilarity(cellchat, type = 'functional')

############## Merge PBMC vs. Tumor ##########

cco.list <- list( PBMC = cco.pbmc,Tumor = cco.tumor)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list),
                          cell.prefix = TRUE)
par(mfrow = c(1,1), xpd = TRUE)
netVisual_diffInteraction(cellchat,
                          weight.scale = T,
                          measure = 'weight',
                          title.name = 'Differential interaction strength between PBMC v.s.T')

netVisual_heatmap(cellchat, measure = 'weight',
                  title.name = 'interaction heatmap between PBMC v.s.T')
num.link <- sapply(cco.list, function(x){
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) #control the dot size in the different datasets
gg <- list()
for(i in 1:length(cco.list)){
  gg[[i]] <- netAnalysis_signalingRole_scatter(cco.list[[i]],
                                               title = names(cco.list)[i],
                                               weight.MinMax = weight.MinMax)}
patchwork::wrap_plots(plots = gg)
pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]],
                                        pattern = 'all',
                                        signaling = pathway.union,
                                        title = names(cco.list)[1],
                                        width = 7.5,height = 14)


ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]],
                                        pattern = 'all',
                                        signaling = pathway.union,
                                        title = names(cco.list)[2],
                                        width = 7.5,height = 14)


draw(ht1 + ht2, ht_gap = unit(1.0, 'cm'))
######################################################
cco.list <- list( Normal = cco.normal,Tumor = cco.tumor)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list),
                          cell.prefix = TRUE)
gg1 <- rankNet(cellchat, mode = 'comparison', stacked = T, do.stat = TRUE)

cco.list <- list( PBMC = cco.pbmc,Tumor = cco.tumor)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list),
                          cell.prefix = TRUE)
gg2 <- rankNet(cellchat, mode = 'comparison', stacked = T, do.stat = TRUE)

gg1+gg2


######################################################
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list , attribute = c('idents','count'))
for(i in 1:length(cco.list)){
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge = F,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions -", names(cco.list)[i]))
}


par(mfrow = c(1,2))
s.cell <- c('CD8+T_1','CD8+T_2', 'B')
count1 <- cco.list[[1]]@net$count[s.cell, s.cell]
count2 <- cco.list[[2]]@net$count[s.cell, s.cell]
weight.max <- max(max(count1), max(count2))
netVisual_circle(count1, weight.scale = T, label.edge = T,
                 edge.weight.max = weight.max, edge.width.max = 12,
                 title.name = paste0('Number of interactions-',names(cco.list)[1]))
netVisual_circle(count2, weight.scale = T, label.edge = T,
                 edge.weight.max = weight.max, edge.width.max = 12,
                 title.name = paste0('Number of intearctions-', names(cco.list)[2]))

gg1 <- rankNet(cellchat, mode = 'comparison', stacked = T, do.stat = TRUE)
gg1


library(ComplexHeatmap)
pathway.union <- union(cco.list[[1]]@netP$pathways,
                       cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]],
                                        pattern = 'all',
                                        signaling = pathway.union,
                                        title = names(cco.list)[1],
                                        width = 10, height = 13)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]],
                                        pattern = 'all',
                                        signaling = pathway.union,
                                        title = names(cco.list)[2],
                                        width = 10, height = 13)
draw(ht1 + ht2, ht_gap = unit(1.0, 'cm'))
#save as compare_signal_pattern_all.pdf 10*6

pathway.union <- union(cco.list[[1]]@netP$pathways,
                       cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]],
                                        pattern = 'outgoing',
                                        signaling = pathway.union,
                                        title = names(cco.list)[1],
                                        width = 10, height = 13)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]],
                                        pattern = 'outgoing',
                                        signaling = pathway.union,
                                        title = names(cco.list)[2],
                                        width = 10, height = 13)
draw(ht1 + ht2, ht_gap = unit(0.5, 'cm'))
#save as compare_signal_pattern_outgoing.pdf


pathway.union <- union(cco.list[[1]]@netP$pathways,
                       cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]],
                                        pattern = 'incoming',
                                        signaling = pathway.union,
                                        title = names(cco.list)[1],
                                        width = 10, height = 13)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]],
                                        pattern = 'incoming',
                                        signaling = pathway.union,
                                        title = names(cco.list)[2],
                                        width = 10, height = 13)
draw(ht1 + ht2, ht_gap = unit(0.5, 'cm'))
#save as compare_signal_pattern_incoming.pdf


pathway.show <- 'MIF'
pathway.show <-'MHC-I'
pathway.show <- 'BTLA'
pathway.show <- 'CD99'
pathway.show <- 'MHC-II'
pathway.show <- 'SELPLG'
pathway.show <- 'ITGB2'

weight.max <- getMaxWeight(cco.list, slot.name = c('netP'),
                           attribute = pathway.show)
par(mfrow = c(1,3), xpd = TRUE)
for (i in 1:length(cco.list)) {
  netVisual_aggregate(cco.list[[i]], signaling = pathway.show,
                      layout = 'circle',
                      edge.weight.max = weight.max[1],
                      edge.width.max = 10,
                      signaling.name = paste0(pathway.show ,'  ', names(cco.list)[i]))
}

#par(mfrow = c(1,2), xpd = TRUE)
#for(i in 1:length(cco.list)){
#  netVisual_aggregate(cco.list[[i]], signaling = pathway.show,
#                      layout = 'chord', pt.title =3,
#                      title.space = 0.05,
#                      vertex.label.cex = 0.5,
#                      signaling.name = paste(pathway.show, names(cco.list)[i]))
#}

levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(4,5),
                      targets.use = c(1,2,3,6),
                      comparison = c(1,2),
                      font.size = 13,
                      font.size.title = 13,
                      angle.x =45)
p

p1 <- netVisual_bubble(cellchat, sources.use = c(4,5),
                       targets.use = c(1,2,3,6),
                       comparison = c(1,2),
                       max.dataset = 2,
                       title.name = 'Increased signaling in TIL',
                       angle.x = 45,
                       remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(4,5),
                       targets.use = c(1,2,3,6), comparison = c(1,2),
                       max.dataset = 1,
                       title.name = 'Decreased signaling in TIL',
                       angle.x = 45,
                       remove.isolate = T)
p <- p1+p2


##compute the network centrality score
cco.tumor <- netAnalysis_computeCentrality(cco.tumor,
                                           slot.name = 'netP')
netAnalysis_signalingRole_network(cco.tumor, signaling = pathway.show,
                                  width = 9, height = 3.5,
                                  font.size = 10)


library(NMF)
selectK(cco.pbmc, pattern = 'outgoing')
nPatterns = 2
cco.pbmc <- identifyCommunicationPatterns(cco.pbmc,
                                          pattern = 'outgoing',
                                          k = nPatterns,
                                          height = 10,
                                          width = 9)
netAnalysis_river(cco.pbmc, pattern = 'outgoing')

selectK(cco.pbmc, pattern = 'incoming')
nPatterns = 2
cco.pbmc <- identifyCommunicationPatterns(cco.pbmc,
                                          pattern = 'incoming',
                                          k = nPatterns,
                                          height = 10,
                                          width = 9)
netAnalysis_river(cco.pbmc, pattern = 'incoming')
###############

selectK(cco.normal, pattern = 'outgoing')
nPatterns = 2
cco.normal <- identifyCommunicationPatterns(cco.normal,
                                            pattern = 'outgoing',
                                            k = nPatterns,
                                            height = 10,
                                            width = 9)
netAnalysis_river(cco.normal, pattern = 'outgoing')

selectK(cco.normal, pattern = 'incoming')
nPatterns = 2
cco.normal <- identifyCommunicationPatterns(cco.normal,
                                            pattern = 'incoming',
                                            k = nPatterns,
                                            height = 10,
                                            width = 9)
netAnalysis_river(cco.normal, pattern = 'incoming')
############################################

selectK(cco.tumor, pattern = 'outgoing')
nPatterns = 2
cco.tumor <- identifyCommunicationPatterns(cco.tumor,
                                           pattern = 'outgoing',
                                           k = nPatterns,
                                           height = 10,
                                           width = 9)
netAnalysis_river(cco.tumor, pattern = 'outgoing')

selectK(cco.tumor, pattern = 'incoming')
nPatterns = 2
cco.tumor <- identifyCommunicationPatterns(cco.tumor,
                                           pattern = 'incoming',
                                           k = nPatterns,
                                           height = 10,
                                           width = 9)
netAnalysis_river(cco.tumor, pattern = 'incoming')


###############
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c('PBMC','Normal','Tumor'))
plotGeneExpression(cellchat, signaling = 'MIF',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'CD99',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'MHC-II',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'SELPLG',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'PECAM1',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'ITGB2',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'ICAM',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'LCK',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'BAG',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'BTLA',
                   split.by = 'datasets',
                   colors.ggplot = T)
plotGeneExpression(cellchat, signaling = 'CD86',
                   split.by = 'datasets',
                   colors.ggplot = T)
