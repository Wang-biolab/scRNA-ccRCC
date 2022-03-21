library(Seurat)
library(dplyr)
library(tidyverse)
library(harmony)
library(patchwork)

pc.num = 1:40
scRNAsub <- RunTSNE(scRNAsub, reduction="harmony", dims=pc.num) %>%
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>%
  FindClusters(resolution=0.1) %>% FindClusters(resolution=0.2) %>%
  FindClusters(resolution=0.3) %>% FindClusters(resolution=0.5)
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "RNA_snn_res.0.1",
              label = T) + ggtitle("RNA_snn_res.0.1")
p2 <- DimPlot(scRNAsub, reduction = "umap", group.by = "RNA_snn_res.0.2",
              label = T) + ggtitle("RNA_snn_res.0.2")
p3 <- DimPlot(scRNAsub, reduction = "umap", group.by = "RNA_snn_res.0.3",
              label = T) + ggtitle("RNA_snn_res.0.3")
p4 <- DimPlot(scRNAsub, reduction = "umap", group.by = "RNA_snn_res.0.5",
              label = T) + ggtitle("RNA_snn_res.0.5")
p <- (p1|p2)/(p3|p4)

scRNAsub$seurat_clusters <- scRNAsub$RNA_snn_res.0.3
Idents(scRNAsub) <- "RNA_snn_res.0.3"

p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "orig")
pc = p1 +  plot_layout(guides = "collect")
save(scRNAsub, file = "Subset_CD8+T3/scRNAsub_harmony.Rdata")
ClusterMarker <- FindAllMarkers(scRNAsub, assay = "RNA", slot = "data", only.pos = T,
                                logfc.threshold = 0.25, min.pct = 0.2)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'Subset_CD8+T3/ClusterMarker.csv', row.names=F)
top = 5   #可根据需要调整
TopMarkers <- ClusterMarker%>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)


markers <- TopMarkers$gene
p <- DoHeatmap(scRNAsub, features = markers, size = 4, label = F)
ggsave("Subset_CD8+T3/Markers_heatmap.png", width = 7, height = 3.5)

p <- VlnPlot(scRNAsub, features = markers, stack = T, split.by = 'orig')
ggsave('Subset_CD8+T3/vlnplot.png',p, width = 7, height = 3.5 )
save(scRNAsub, file = "Subset_CD8+T3/scRNAsub_harmony.Rdata")
