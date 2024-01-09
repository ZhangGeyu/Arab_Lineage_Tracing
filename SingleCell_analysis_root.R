library('Seurat')

lsSeurat <- list()
for (i in c('seedling2_root')){
  path = paste0("/",i,"/output/filter_matrix/")
  raw.data <- Read10X(data.dir = path, gene.column = 1)
  raw.data <- raw.data[grep('Zea|4sgRNA',rownames(raw.data)),]
  rownames(raw.data) <- gsub("Plant__ZeawithTE_","",rownames(raw.data))
  raw.data[1:5,1:5]
  dim(raw.data)
  dim(raw.data)
  my.seurat <- CreateSeuratObject(counts = raw.data)
  my.seurat <- NormalizeData(my.seurat, normalization.method = "LogNormalize",
                             scale.factor = 1e6)
  my.seurat <- SCTransform(my.seurat)
  my.seurat <- RunPCA(object = my.seurat, verbose = FALSE, npcs = 100)
  my.seurat <- RunUMAP(object = my.seurat, dims = 1:30, 
                       min.dist = 0.05, n.neighbors = 5, seed.use = 100)
  lsSeurat[[i]] <- my.seurat
}
names(lsSeurat)
gene <- gsub('_','-',rownames(raw.data)[grep('4sgRNA',rownames(raw.data))])

for (i in names(lsSeurat)){
  print(i)
  DefaultAssay(lsSeurat[[i]]) <- 'RNA'
  print(FeaturePlot(lsSeurat[[i]],features = gene,ncol = 5))
  DefaultAssay(lsSeurat[[i]]) <- 'SCT'
}

for (i in names(lsSeurat)){
  print(dim(lsSeurat[[i]]@assays$RNA@counts))
}

datCellId <- data.frame()
for (i in names(lsSeurat)){
  tmp <- data.frame(cell=colnames(lsSeurat[[i]]@assays$RNA@counts))
  tmp$sample <- i
  datCellId <- rbind(datCellId,tmp)
}
table(datCellId$sample)
write.table(datCellId,file = 'cellID-root.tsv',sep = '\t',row.names = F,quote = F)
getwd()

maize.integrated <- lsSeurat$`seedling2_root`
maize.integrated <- RunPCA(object = maize.integrated, verbose = FALSE, npcs = 100)
maize.integrated <- RunUMAP(object = maize.integrated, dims = 1:30)#, min.dist = 0.05, n.neighbors = 5, seed.use = 100)
head(maize.integrated@meta.data)

maize.integrated@meta.data$sampleid <- 1

maize.integrated@meta.data$sample <- 'N'
for (i in 1:length(names(lsSeurat))){
  maize.integrated@meta.data[maize.integrated@meta.data$sampleid==i,]$sample <- names(lsSeurat)[i]
}
table(maize.integrated@meta.data$sample)
DimPlot(maize.integrated, reduction = "umap",split.by = 'sample',group.by = 'sample',)


maize.integrated <- FindNeighbors(maize.integrated, reduction = "pca", dims = 1:30)
maize.integrated <- FindClusters(maize.integrated, resolution = 0.5)
DimPlot(maize.integrated, reduction = "umap",label = T)

# Plasmid express level

DefaultAssay(maize.integrated) <- "RNA"
maize.integrated <- NormalizeData(maize.integrated, normalization.method = "LogNormalize",
                                  scale.factor = 1e6)
gene <- gsub('_','-',rownames(raw.data)[grep('4sgRNA',rownames(raw.data))])
png(filename = paste0("Plasmid_express.png"),width = 2000,height = 8000,res = 200)
FeaturePlot(maize.integrated,features = gene,split.by = 'sample',ncol = 4)
dev.off()

# TE express level

TE <- gsub('_','-',c('TE_homo_1053200','TE_homo_1053202','TE_homo_1053284','TE_homo_1053359',
             'TE_homo_1053283','TE_homo_1053355','TE_homo_1053358','TE_homo_1053188'))
png(filename = paste0("TE_express.png"), width = 2000,height = 6000,res = 200)
FeaturePlot(maize.integrated,features = TE,split.by = 'sample',ncol = 4)
dev.off()


# cell type annotation

library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)

gene <- as.data.frame(fread("",header = TRUE)) # File of marker genes
head(gene)
gene_list <- gene$gene

table(gene$cell_type)
gene <- gene %>% group_by(gene) %>% summarise(celltype=paste(cell_type,collapse = ", ")) %>% as.data.frame()
gene$name <- paste0(gene$celltype,"-",gene$gene)
rownames(gene) <- gene$gene


my.seurat <- maize.integrated
head(my.seurat@meta.data)
str(my.seurat@meta.data)
Idents(my.seurat) <- "seurat_clusters"
table(Idents(my.seurat))
table(my.seurat@meta.data$seurat_clusters)
head(gene)

selected_gene_indices <- which(rownames(my.seurat@assays$RNA@data) %in% gene$gene)
my.seurat@assays$RNA@data <- my.seurat@assays$RNA@data[selected_gene_indices, ]


my.seurat@assays$RNA@data[1:5,1:5]
rownames(my.seurat@assays$RNA@data) <- gene[rownames(my.seurat@assays$RNA@data),]$name
DefaultAssay(my.seurat) <- "RNA"

features <- gene$name
DotPlot(my.seurat, features = features) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))


plt.maize <- as.data.frame(maize.integrated@reductions$umap@cell.embeddings)
head(plt.maize)
colnames(plt.maize) <- c('V1','V2')
plt.maize <- merge(plt.maize,maize.integrated@meta.data,by='row.names')
rownames(plt.maize) <- plt.maize$Row.names
plt.maize$cell <- rownames(plt.maize)
plt.maize$cell <- paste0(plt.maize$cell,'_',plt.maize$sample)
head(plt.maize)

write.table(plt.maize,'plt_maize_root.txt', sep = '\t', quote = FALSE)

