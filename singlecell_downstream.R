directory = '/path/to/expression/matrix/'
setwd(directory)

# ggplot2 v3.4.2
library(ggplot2)
# dplyr v1.1.2
library(dplyr)
# Seurat v4.3.0.1
library(Seurat)

scdata <- Read10X(data.dir = './data')
scdata <- CreateSeuratObject(counts = scdata, project = 'G4TotiSC')

tmp <- scdata
scdata <- tmp

scdata[['percent.mt']] <- PercentageFeatureSet(scdata, pattern = '^mt-')
scdata <- subset(scdata, subset = nFeature_RNA >= 1000 & percent.mt <= 10)

scdata <- NormalizeData(scdata, normalization.method = 'LogNormalize', scale.factor = 10000)
scdata <- FindVariableFeatures(scdata, selection.method = 'vst', nfeatures = 2000)
scdata <- ScaleData(scdata, features = row.names(scdata))

scdata <- RunPCA(scdata, features = VariableFeatures(object = scdata))
scdata <- FindNeighbors(scdata)
scdata <- FindClusters(scdata, resolution = 0.63)
scdata.markers <- FindAllMarkers(scdata, only.pos = TRUE, min.pct = 0.25)
write.csv(scdata.markers, file = './18cluster_marker.csv', na = ' ')

# scdata <- RunTSNE(scdata, dims = 1:10)
# DimPlot(scdata, reduction = 'tsne', label = TRUE) + NoLegend()
scdata <- RunUMAP(scdata, dims = 1:10)
DimPlot(scdata, reduction = 'umap', label = TRUE) + NoLegend()

topn <- scdata.markers %>% group_by(cluster) %>% top_n(n = 6, wt = avg_log2FC)
DoHeatmap(scdata, features = topn$gene)

result <- scdata
result.markers <- scdata.markers


number <- c(8, 17, 5, 13, 9, 12, 7, 11, 2, 6, 14, 15, 1, 16, 18, 3, 10, 4)
new.cluster.ids <- number
names(new.cluster.ids) <- levels(scdata)
scdata <- RenameIdents(scdata, new.cluster.ids)
DimPlot(scdata, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()
ggsave('umap.svg', width = 8, height = 8, dpi = 300)

order <- c(
  '1', '2', '3', '4', '5', '6', '7', '8', '9', 
  '10', '11', '12', '13', '14', '15', '16', '17', '18'
)
levels(scdata) <- levels(factor(x = levels(scdata), levels = order))
gene <- read.csv('gene_heatmap.csv', header = FALSE, col.names = c("gene"))
DoHeatmap(scdata, features = gene$gene, label = FALSE) + scale_fill_gradientn(colours = c("#286cae", "white", "#b6223d"))
ggsave("heatmap.svg", width = 8, height = 7, dpi = 300)

DoHeatmap(AverageExpression(scdata, return.seurat = TRUE), features = gene$gene, label = FALSE, draw.lines = FALSE) + scale_fill_gradientn(colours = c("#286cae", "white", "#b6223d"))
ggsave("heatmap_avg.svg", width = 4, height = 7, dpi = 300)

DotPlot(scdata, features = gene$gene, assay = 'RNA') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
ggsave("dotplot.svg", width = 16, height = 6, dpi = 300)

expr_level <- AverageExpression(subset(scdata, features = gene$gene))
write.csv(expr_level, 'expression_level.csv')


# tidyr v1.3.0
library(tidyr)
# paletteer v1.6.0
library(paletteer)

target <- read.csv('gene_violin.csv', header = FALSE, col.names = c("gene"))
vln.df <- subset(scdata, features = target$gene)@assays$RNA@scale.data %>%
  t() %>%
  as.data.frame()%>%
  tibble::rownames_to_column("CB") %>% 
  mutate(cluster = factor(number[scdata$seurat_clusters]))%>%
  pivot_longer(cols = 2:(ncol(.)-1), names_to = "gene", values_to = "exp")%>% 
  mutate(gene = factor(gene, levels = target$gene))

vlcolor = paletteer_d(`"ggsci::default_nejm"`)
color = colorRampPalette(vlcolor)(length(unique(vln.df$cluster)))

ggplot(vln.df, aes(cluster, exp), color = cluster) +
  geom_violin(aes(fill = cluster), scale = "width") +
  scale_fill_manual(values = color) +
  facet_grid(gene~., scales = "free_y") +
  scale_y_continuous(expand = c(0, 0.5)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom = element_text(hjust = NULL, vjust = NULL, color = "black", size = 18),
    axis.title.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    legend.position = "none",
    panel.spacing.y = unit(0, "cm"),
    strip.text.y = element_text(angle=0, size = 18, hjust = 0),
    strip.background.y = element_blank()
  )
ggsave('violin.svg', width = 10, height = 10, dpi = 300)
