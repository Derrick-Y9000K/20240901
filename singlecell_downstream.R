# R_v4.2.2
# single cell downstream analysis

# Seurat_v4.3.0.1
# SeuratObject_v4.1.3
library(Seurat)
# ggplot2_v3.4.2
library(ggplot2)
# dplyr_v1.1.2
library(dplyr)
# tidyr_v1.3.0
library(tidyr)
# paletteer_v1.6.0
library(paletteer)

directory = '/path/to/PlaYok/directory'
setwd(directory)

scdata <- Read10X(data.dir = './data')
scdata <- CreateSeuratObject(counts = scdata, project = 'PlaYok')

scdata[['percent.mt']] <- PercentageFeatureSet(scdata, pattern = '^mt-')
scdata <- subset(scdata, subset = nFeature_RNA >= 1000 & percent.mt <= 10)

scdata <- NormalizeData(scdata, normalization.method = 'LogNormalize', scale.factor = 10000)
scdata <- FindVariableFeatures(scdata, selection.method = 'vst', nfeatures = 2000)
scdata <- ScaleData(scdata)

scdata <- RunPCA(scdata, features = VariableFeatures(object = scdata))
scdata <- FindNeighbors(scdata)
scdata <- FindClusters(scdata, resolution = 0.63)
scdata <- RunUMAP(scdata, dims = 1:10)
DimPlot(scdata, reduction = 'umap', label = TRUE, pt.size = 0.5)

scdata.markers <- FindAllMarkers(scdata, only.pos = TRUE, min.pct = 0.25)
write.csv(scdata.markers, file = './marker.csv', na = ' ')

number <- c(7, 17, 11, 13, 8, 12, 6, 10, 2, 5, 14, 15, 1, 16, 18, 3, 9, 4)
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
target <- read.csv('expression_level_norm-0613.csv', header = TRUE)$gene
DoHeatmap(scdata, features = target, label = FALSE) + scale_fill_gradientn(colours = c('#286cae', 'white', '#b6223d'))
ggsave('heatmap.svg', width = 8, height = 7, dpi = 300)

DotPlot(scdata, features = target, assay = 'RNA') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
ggsave('dotplot.svg', width = 16, height = 6, dpi = 300)

target <- read.csv('genes.csv', header = FALSE, col.names = c('gene'))$gene
vln.df <- scdata@assays$RNA@scale.data %>%
  t() %>%
  as.data.frame()%>%
  select(target) %>% 
  tibble::rownames_to_column('CB') %>% 
  mutate(cluster = factor(number[scdata$seurat_clusters]))%>%
  pivot_longer(cols = 2:(ncol(.)-1), names_to = 'gene', values_to = 'exp')%>% 
  mutate(gene = factor(gene, levels = target))

vlcolor = paletteer_d(`"ggsci::default_nejm"`)
vlcolor = colorRampPalette(vlcolor)(length(unique(vln.df$cluster)))

ggplot(vln.df, aes(cluster, exp), color = cluster) +
  geom_violin(aes(fill = cluster), scale = 'width') +
  scale_fill_manual(values = vlcolor) +
  facet_grid(gene~., scales = 'free_y') +
  scale_y_continuous(expand = c(0, 0.5)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title.x.bottom = element_blank(),
    axis.ticks.x.bottom = element_blank(),
    axis.text.x.bottom = element_text(hjust = NULL, vjust = NULL, color = 'black', size = 18),
    axis.title.y.left = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left = element_blank(),
    legend.position = 'none',
    panel.spacing.y = unit(0, 'cm'),
    strip.text.y = element_text(angle=0, size = 18, hjust = 0),
    strip.background.y = element_blank()
  )
ggsave('violin.svg', width = 10, height = 10, dpi = 300)
