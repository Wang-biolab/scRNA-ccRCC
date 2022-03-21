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
load("Subset_CD8+T2/scRNAsub_harmony.Rdata")

geneList = deg1$avg_log2FC
names(geneList) = deg1$gene
geneList = sort(geneList, decreasing = T)
geneList[1:10]

geneset <- read.gmt('h.all.v7.4.symbols.gmt')
geneset$term = str_remove(geneset$term, 'HALLMARK_')
head(geneset)


expr <- as.data.frame(scRNAsub@assays$RNA@data)
expr[1:3,1:4]
expr = as.matrix(expr)

deg <- ClusterMarker
item < intersect(rownames(deg), rownames(expr))
expr  <- expr[item,]

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


deg <- FindMarkers(scRNAsub, ident.1 = c('0','1','4','5'), ident.2 = c('2','3'))

#deg <- FindMarkers(scRNAsub, ident.1 = 'T', ident.2 = c('N','P'))
#deg <- FindAllMarkers(scRNAsub, assay = "RNA", slot = "data", only.pos = T,
  #                           logfc.threshold = 0.25, min.pct = 0.2)
deg1 <- data.frame(gene = rownames(deg), deg)
head(deg)

geneList = deg1$avg_log2FC
names(geneList) = deg1$gene
geneList = sort(geneList, decreasing = T)
geneList[1:10]

egmt <- GSEA(geneList, TERM2GENE = geneset, verbose = F, pvalueCutoff = 0.1)
y=data.frame(egmt)
head(y)
write.csv(y, file = 'Subset_CD8+T2/GSEA_tumor.enrichment.csv', quote = F)
dotplot(egmt, split = '.sign')+facet_grid(~.sign)
gseaplot2(egmt, geneSetID = 1, title = egmt$Description[1])








log_FC_t = 0.5
P.Value_t = 0.05
k1 = (deg1$p_val_adj < P.Value_t)&(deg1$avg_log2FC < -log_FC_t)
k2 = (deg1$p_val_adj < P.Value_t)&(deg1$avg_log2FC > log_FC_t)
table(k1)
table(k2)

change = ifelse(k1, 'down', ifelse(k2, 'up', 'stable'))
deg1$change <- change
head(deg1)

s2e <- bitr(deg1$gene,
            fromType = 'SYMBOL',
            toType = 'ENTREZID',
            OrgDb = org.Hs.eg.db)
deg1 <- inner_join(deg1, s2e, by= c('gene'='SYMBOL'))
head(deg1)

gene_up = deg1[deg1$change == 'up','gene']
gene_down = deg1[deg1$change == 'down', 'gene']
gene_diff = c(gene_up, gene_down)

gene_all = deg1[,'ENTREZID']
gene_up_KEGG = deg1[deg1$change == 'up','ENTREZID']
gene_down_KEGG = deg1[deg1$chang == 'down', 'ENTREZID']
gene_diff_KEGG = c(gene_up_KEGG, gene_down_KEGG)

ego_BP <- enrichGO(gene = gene_up,
                   keyType = 'SYMBOL',
                   OrgDb = org.Hs.eg.db,
                   ont = 'CC',
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
head(summary(ego_BP))
dotplot(ego_BP, showCategory = 10)

kk.up <- enrichKEGG(gene = gene_up_KEGG,
                    organism = 'hsa',
                    universe = gene_all,
                    pvalueCutoff = 0.5,
                    qvalueCutoff = 0.5
                    )
kk.down <- enrichKEGG(gene = gene_down_KEGG,
                      organism = 'hsa',
                      universe = gene_all,
                      pvalueCutoff = 0.5,
                      qvalueCutoff = 0.5)
ekegg <-setReadable(kk.up, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
kegg_diff_dt <- data.frame(ekegg)
head(kegg_diff_dt)

dotplot(ekegg, showCategory = 10)
write.csv(kegg_diff_dt,file = 'Subset_CD8+T2/kegg_diff_in_CD8+T2_tumor.csv')
#####################