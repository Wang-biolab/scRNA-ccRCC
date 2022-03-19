############## Monocle analysis #####################3
rm(list = ls())
load("Subset_CD8+T3/scRNAsub_harmony.Rdata")
Idents(scRNAsub) <- 'orig'
scRNAsub.T <- subset(scRNAsub, idents = 'T')
dir.create('Subset_CD8+T3/Monocle')
library(monocle)

data <- GetAssayData(scRNAsub.T, assay = 'RNA', slot = 'counts')
pd <- new('AnnotatedDataFrame', data = scRNAsub.T@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds)

order.genes <- SCTransform(scRNAsub.T) %>% VariableFeatures()
mycds <- setOrderingFilter(mycds, order.genes)

mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds)

plot1 <- plot_cell_trajectory(mycds, color_by = 'State')
ggsave('Subset_CD8+T3/CD8_T3_trajectory_state.png', plot1, width = 5, height = 3.5)

plot2 <- plot_cell_trajectory(mycds, color_by = 'seurat_clusters')
ggsave('Subset_CD8+T3/CD8_T3_trajectory_seuratClusters.png',plot2, width = 5, height = 3.5)


plot3 <- plot_cell_trajectory(mycds, color_by = 'Pseudotime')
ggsave('Subset_CD8+T3/CD8_T3_trajectory_pseudotime.png', plot3, width = 5, height = 3.5)

plotc <- plot2|plot1|plot3
ggsave('Subset_CD8+T3/CD8_T3_trajectory_combination.png', plotc, width = 15.4, height = 4.5)

save(mycds, file = 'Subset_CD8+T3/CD8_T3.mycds.Rdata')

s.genes <- c('LAG3','GZMK','FXYD2','TRBV20-1','TRDV1')
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = 'State', color_by = 'State')
p2 <- plot_genes_violin(mycds[s.genes,], grouping = 'State', color_by = 'State')
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = 'State')
ggsave(file = 'Subset_CD8+T3/state3_genes_jitterplot.png', p1, width = 7.5, height = 5)
ggsave(file = 'Subset_CD8+T3/state3_genes_violinplot.png', p2, width = 7.5, height = 5)
ggsave(file = 'Subset_CD8+T3/state3_genes_pseudotime.png', p3, width = 5, height = 4)

s.genes <- c('GNLY','FGFBP2','SPON2','GZMB','RPS28')
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = 'State')
ggsave(file = 'Subset_CD8+T3/genes_pseudotime234.png', p3, width = 5, height = 4)
s.genes <- c('LYZ','CST3','IL1B','GPX1','C1QC')
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = 'State')
ggsave(file = 'Subset_CD8+T3/genes_pseudotime1.png', p3, width = 5, height = 4)

Time_diff <- differentialGeneTest(mycds[order.genes,], cores = 1,
                                  fullModelFormulaStr = '~sm.ns(Pseudotime)')

Time_diff <- Time_diff[,c(5,2,3,4,1,6)]
write.csv(Time_diff, file = 'Subset_CD8+T3/Time_diff.csv', row.names = F)

Time_genes <- Time_diff[order(Time_diff$qval, decreasing = F)[1:50],1]
p = plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters = 3,
                            show_rownames = T, return_heatmap = T)
ggsave('Subset_CD8+T3/Time_heatmap.png', p, width = 8.5, height = 6)

hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c('gene_short_name', 'pval','qval')]

beam_res <- BEAM(mycds[order.genes,], branch_point = 1)
beam_res <- beam_res[, c(5,2,3,4,1,6)]
BEAM_genes <- beam_res[order(beam_res$qval, decreasing = F)[1:50],1]

p <- plot_genes_branched_heatmap(mycds[BEAM_genes, ], branch_point = 1,
                                 num_clusters = 2, show_rownames = T, return_heatmap = T)
ggsave('Subset_CD8+T3/BEAM_heatmap.png', p$ph_res, width = 8.5, height = 6)

beam_res <- BEAM(mycds[order.genes,], branch_point = 2)
beam_res <- beam_res[, c(5,2,3,4,1,6)]
BEAM_genes <- beam_res[order(beam_res$qval, decreasing = F)[1:50],1]
p <- plot_genes_branched_heatmap(mycds[BEAM_genes, ], branch_point = 2,
                                 num_clusters = 2, show_rownames = T, return_heatmap = T)
ggsave('Subset_CD8+T3/BEAM_heatmap_branch2.png', p$ph_res, width = 8.5, height = 6)

hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- beam_res[hp.genes, c('gene_short_name', 'pval','qval')]
write.csv(BEAM_sig, 'Subset_CD8+T3/BEAM_sig.csv', row.names = F, quote = F)
