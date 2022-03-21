###GSEA & GSVA enrichment analysis ##
library(Seurat)
library(GSVA)
library(dplyr)
library(GSEABase)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(enrichplot)
library(clusterProfiler)
options(stringAsFactor = F)
load("Subset_CD8+T3/scRNAsub_harmony.Rdata")
geneset <- read.gmt('h.all.v7.4.symbols.gmt')
geneset$term = str_remove(geneset$term, 'HALLMARK_')
head(geneset)

deg <- FindAllMarkers(scRNAsub, assay = "RNA", slot = "data", only.pos = T,
                      logfc.threshold = 0.25, min.pct = 0.2)

expr <- as.data.frame(scRNAsub@assays$RNA@data)
expr[1:3,1:4]
expr = as.matrix(expr)
item <- intersect(rownames(deg), rownames(expr))
expr <- expr[item, ]
kegg_list = split(geneset$gene, geneset$term)
kegg_list[1:2]

kegg1 <- gsva(expr, kegg_list,
              kcdf = 'Gaussian',
              method = 'gsva',
              parallel.sz = 12)

kegg2 <- kegg1
kegg2[1:10,1:4]
meta <- as.data.frame(scRNAsub@meta.data[, c('orig','seurat_clusters')])

meta <- meta %>% arrange(meta$seurat_clusters)
kegg2 <- kegg2[, rownames(meta)]
identical(colnames(kegg2), rownames(meta))

kegg3 <- cbind(meta, t(kegg2)) %>% tibble::rownames_to_column()
kegg4 <- tidyr::pivot_longer(kegg3, cols = 4:ncol(kegg3), 'kegg', 'value')
kegg5 <- dplyr::group_by(kegg4, seurat_clusters, kegg) %>% dplyr::summarise(value2 = mean(value))
kegg6 <- tidyr::pivot_wider(kegg5, names_from = 'kegg', values_from = 'value2') %>%
  tibble::column_to_rownames(var = 'seurat_clusters')

library(pheatmap)
G1 <- pheatmap(kegg6,
               cluster_rows = T,
               cluster_cols = T,
               show_rownames = T,
               show_colnames = T,
               color = colorRampPalette(c('blue','white','red'))(100),
               cellwidth = 10,
               cellheight = 14,
               fontsize = 10)


#markers <- c('LAG3','CTLA4','FXYD2','TRBV20-1','PDCD1'
#             )

#p <- DoHeatmap(scRNAsub, features = markers, size = 4, label = F)
ggsave("Subset_CD8+T3/Markers_heatmap.png", width = 7, height = 3.5)

