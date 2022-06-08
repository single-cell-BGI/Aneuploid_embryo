library(dplyr)
library(Seurat)
library(monocle)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library("viridis")
library(reshape2)

MTN.color = brewer.pal(8,"Dark2")[1:3]
names(MTN.color) = c("E","T","M")
load("../CNVfiltered_11561cells_combine_BGI500DipseqSeed30res017.updatemeta.RData")
EPI = subset(embryo.integrated_CNVdefined,subset=ChangeNNcelltype %in% c("EPI"))
chrom.state = c("E",paste0(c(1:22,"X","Y"),"T"),paste0(c(1:22,"X"),"M"))

EPI <- CreateSeuratObject(counts = as.matrix(EPI@assays$RNA@counts),meta.data = EPI@meta.data,project = "EPI", min.cells = 1, min.features = 200)
EPI <- NormalizeData(EPI)
EPI <- FindVariableFeatures(EPI, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(EPI)

#######cell cycle######
get.id = function(y){y = paste0("^",y,"-ENSG");grep(pattern = y,x = rownames(EPI@assays$RNA@counts),value = T)}
get.id = function(y){y = paste0("^",y,"-ENSG");grep(pattern = y,x = rownames(embryo.integrated_CNVdefined@assays$RNA@counts),value = T)}

s.genes <- cc.genes$s.genes
s.genes = unlist(lapply(as.list(s.genes),get.id))
g2m.genes <- cc.genes$g2m.genes
g2m.genes = unlist(lapply(as.list(g2m.genes),get.id))
EPI <- CellCycleScoring(EPI, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)

EPI <- ScaleData(EPI,features = all.genes) #, vars.to.regress = c("nCount_RNA", "percent.mt")

EPI <- RunPCA(EPI, features = VariableFeatures(object = EPI))
EPI <- JackStraw(EPI, num.replicate = 100,dims = 40)
dim = 12 ###18
EPI <- ScoreJackStraw(EPI, dims = 1:dim)
ElbowPlot(EPI,ndims = 40)
EPI <- FindNeighbors(EPI, dims = 1:dim)
EPI<- FindClusters(EPI, resolution = 0.3,n.start = 5,n.iter = 1000)

EPI <- RunTSNE(EPI, dims = 1:20,perplexity = 60)
EPI <- RunUMAP(EPI,n.neighbors = 20L, dims = 1:20,min.dist = 0.2,negative.sample.rate = 1)

DimPlot(EPI, label = TRUE,group.by = "seurat_clusters",reduction = "tsne")

EPI@meta.data$seurat_clusters = as.character(EPI@meta.data$seurat_clusters)
EPI@meta.data$seurat_clusters[grep("3",EPI@meta.data$seurat_clusters)] = "C4"
EPI@meta.data$seurat_clusters[grep("2",EPI@meta.data$seurat_clusters)] = "C3"
EPI@meta.data$seurat_clusters[grep("1",EPI@meta.data$seurat_clusters)] = "C2"
EPI@meta.data$seurat_clusters[grep("0",EPI@meta.data$seurat_clusters)] = "C1"
EPI@meta.data$seurat_clusters[grep("C3",EPI@meta.data$seurat_clusters)] = "CT"
EPI@meta.data$seurat_clusters[grep("C1",EPI@meta.data$seurat_clusters)] = "C3"
EPI@meta.data$seurat_clusters[grep("CT",EPI@meta.data$seurat_clusters)] = "C1"

EPI@meta.data$chrom.state = gsub("NN","E",EPI@meta.data$chrom.state)
EPI@meta.data$ab.type = gsub("N","E",EPI@meta.data$ab.type)

########draw cells information in tsne plot
p1 = DimPlot(EPI,group.by = "seurat_clusters",cols = color,reduction = "tsne",pt.size = 1.5)
p2 = DimPlot(EPI,group.by = "ab.type",cols = MTN.color,reduction = "tsne",pt.size = 1.5)
p3 = DimPlot(EPI,group.by = "platform",cols = color[8:9],reduction = "tsne",pt.size = 1.5)
p4 = DimPlot(EPI,group.by = "Phase",cols = color[5:7],reduction = "tsne",pt.size = 1.5)
pdf("EPI_resubcluster.res0.3.dim12.tsne.pdf",width = 8,height = 6)
plot_grid(p1,p2,p3,p4,ncol = 2)
dev.off()
########draw cells information in umap plot
p1 = DimPlot(EPI,group.by = "seurat_clusters",cols = color,reduction = "umap",pt.size = 1.5)
p2 = DimPlot(EPI,group.by = "ab.type",cols = MTN.color,reduction = "umap",pt.size = 1.5)
p3 = DimPlot(EPI,group.by = "platform",cols = color[8:9],reduction = "umap",pt.size = 1.5)
p4 = DimPlot(EPI,group.by = "Phase",cols = color[5:7],reduction = "umap",pt.size = 1.5)

pdf("EPI_resubcluster.res0.3.dim12.umap.pdf",width = 8,height = 6)
plot_grid(p1,p2,p3,p4,ncol = 2)
dev.off()

saveRDS(EPI,"15.EPI.subcluster/EPI.resubcluster.4cluster.rds")

############summary the distribution of aneuploidy cells in subclusters#########
chrom.state = c("E",paste0(c(1:22,"X","Y"),"T"),paste0(c(1:22,"X"),"M"))
EPI@meta.data$chrom.state = factor(EPI@meta.data$chrom.state,levels = chrom.state)
temp = data.frame(table(EPI@meta.data$seurat_clusters,EPI@meta.data$chrom.state),stringsAsFactors = F)
#temp$Var2 = gsub("NN","E",temp$Var2)
colnames(temp) = c("cluster","chrom.state","cell.number")
temp$chrom.state = factor(temp$chrom.state,levels = chrom.state)
tt = table(EPI@meta.data$chrom.state)
tt = tt[tt>2]
temp = temp[which(temp$chrom.state %in% names(tt)),]
temp$cluster = factor(temp$cluster,levels = c("C1","C2","C4","C3"))
p1 = ggplot(data = temp, mapping = aes(x = chrom.state, y = cell.number, fill = cluster)) + 
  geom_bar(stat = 'identity', position = 'fill')+theme_bw()+scale_fill_manual(values = color)+
  theme(axis.text.x = element_text(size=8))+theme_classic()

pdf("EPI_subcluster_chrom.state.dist.rm2cells.pdf",width = 12,height = 6)
print(p1)
dev.off()

color = brewer.pal(9, "Set1")[1:4]
names(color) = c("C1","C2","C3","C4")

temp = data.frame(table(EPI@meta.data$seurat_clusters,EPI@meta.data$ab.type))
colnames(temp) = c("cluster","ab.type","cell.number")
temp$cluster = factor(temp$cluster,levels = c("C1","C2","C4","C3"))
temp$ab.type = factor(temp$ab.type,levels = c("E","T","M"))
p1 = ggplot(data = temp, mapping = aes(x = cluster, y = cell.number, fill = ab.type)) + 
  geom_bar(stat = 'identity', position = 'fill')+theme_bw()+scale_fill_manual(values = MTN.color)+
  theme(axis.text.x = element_text(size=8))+theme_classic()
p2 = ggplot(data = temp, mapping = aes(x = cluster, y = cell.number, fill = ab.type)) + 
  geom_bar(stat = 'identity', position = 'stack')+theme_bw()+scale_fill_manual(values = MTN.color)+
  theme(axis.text.x = element_text(size=8))+theme_classic()

p3 = ggplot(data = temp, mapping = aes(x = ab.type, y = cell.number, fill = cluster)) + 
  geom_bar(stat = 'identity', position = 'fill')+theme_bw()+scale_fill_manual(values = color)+
  theme(axis.text.x = element_text(size=8))+theme_classic()
p4 = ggplot(data = temp, mapping = aes(x = ab.type, y = cell.number, fill = cluster)) + 
  geom_bar(stat = 'identity', position = 'stack')+theme_bw()+scale_fill_manual(values = color)+
  theme(axis.text.x = element_text(size=8))+theme_classic()

pdf("15.EPI.resubcluster/EPI_subcluster_ab.type.dist.pdf",width = 5,height = 4)
plot_grid(p1,p2,p3,p4,ncol = 2)
dev.off()


##########find subcluster-specific genes
Idents(EPI) = "seurat_clusters"
EPI.markers <- FindAllMarkers(EPI, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,return.thresh = 0.05)
EPI.markers.filter = EPI.markers[which(EPI.markers$avg_logFC>0.25 & EPI.markers$pct.1>0.25 & EPI.markers$p_val_adj<0.05),]

##########draw heatmap for subcluster-specific genes
EPI.markers.filter = EPI.markers.filter[order(EPI.markers.filter$cluster,-EPI.markers.filter$avg_logFC),]
DEgenes.mtx = as.matrix(EPI@assays$RNA@counts)[as.matrix(EPI.markers.filter$gene),]
rownames(DEgenes.mtx) = EPI.markers.all.mean.filter$gene.name
meta = EPI@meta.data
meta$ab.type = factor(meta$ab.type,levels = c("E","T","M"))
meta = meta[order(meta$seurat_clusters,meta$ab.type),]
DEgenes.mtx = DEgenes.mtx[,rownames(meta)]

library(ComplexHeatmap)
col1 = colorRampPalette(c("#b2182b","#f4a582"))(12)
col2 = colorRampPalette(c("#fddbc7","#d1e5f0"))(2)
col3 = colorRampPalette(c("#92c5de", "#2166ac"))(12)
colours = rev(c(col1,col2,col3))
mat_scaled = t(scale(t(log1p(DEgenes.mtx)),center = T))
mat_scaled = mat_scaled/max(abs(mat_scaled))*4
mat_scaled[mat_scaled>2] = 2
mat_scaled[mat_scaled< -2] = -2

sub.color = brewer.pal(9, "Set1")[1:4]
names(sub.color) = paste0("C",c(1:4))
ha = HeatmapAnnotation(
  subcluster = c(as.character(meta$seurat_clusters)), 
  ab.type = c(as.character(meta$ab.type)),
  col = list(subcluster = sub.color,
             ab.type = c(MTN.color)),
  gap = unit(c(2, 2), "mm"))
pdf("EPI_subcluster_all_DEgenes.logFC0.25.adjp0.05.complexheatmap.pdf",width = 8,height = 10)
Heatmap(mat_scaled, name = "Z-score", row_split = EPI.markers.all.mean.filter$cluster, column_split = as.character(meta$seurat_clusters),col = colours,show_column_names=F,row_names_gp = gpar(fontsize = 1),cluster_row_slices = F,cluster_column_slices = F,top_annotation = ha,row_dend_width = unit(12, "mm"),row_dend_gp = gpar(lwd = 0.1),show_row_dend = F,show_column_dend = F)
dev.off()

###############GO term#############
library(clusterProfiler)  
library(org.Hs.eg.db)

EPI.markers = EPI.markers.all.mean.filter
combine.count =data.frame(Description = "NA")
combine.p =data.frame(Description = "NA")
combine.gene = data.frame(Description = "NA")
combine.padj = data.frame(Description = "NA")
for(i in paste0("C",c(1:4))){
  temp = EPI.markers[EPI.markers$cluster==i,]
  gene = as.character(temp$gene.name)
  mygene<-select(org.Hs.eg.db,columns=c("SYMBOL","ENTREZID"),keytype="SYMBOL",keys=gene)
  mygene<-mygene$ENTREZID
  ego<-enrichGO(OrgDb = "org.Hs.eg.db",gene=mygene,ont="BP",pvalueCutoff=0.05,readable=TRUE)
  temp = ego@result[ego@result$pvalue<0.05,]
  temp1 = temp[,c("Description","pvalue")]
  temp2 = temp[,c("Description","Count")]
  temp3 = temp[,c("Description","geneID")]
  temp4 = temp[,c("Description","p.adjust")]
  combine.count = merge(combine.count,temp2,by = "Description",all=T)
  combine.p = merge(combine.p,temp1,by = "Description",all=T)
  combine.gene = merge(combine.gene,temp3,by = "Description",all=T)
  combine.padj = merge(combine.padj,temp4,by = "Description",all=T)
}
rownames(combine.p) = combine.p$Description
rownames(combine.padj) = combine.padj$Description
rownames(combine.count) = combine.count$Description
rownames(combine.gene) = combine.gene$Description
combine.gene = combine.gene[-1,-1]
combine.p = combine.p[-1,-1]
combine.count = combine.count[-1,-1]
combine.padj = combine.padj[-1,-1]
colnames(combine.p) = paste0("C",c(1:4))
colnames(combine.count) = paste0("C",c(1:4))
colnames(combine.padj) = paste0("C",c(1:4))
colnames(combine.gene) = paste0("C",c(1:4))


##########EPI trajectory analysis####################
library(monocle)
count_matrix<-EPI@assays$RNA@counts
count_matrix = count_matrix[rowSums(as.matrix(count_matrix))>0,]
sample_sheet<-EPI@meta.data
gene_annotation<-data.frame("gene_short_name"=rownames(count_matrix),"gene_name"=rownames(count_matrix))
row.names(gene_annotation)<-gene_annotation$gene_short_name
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
EPI.trajectory <- newCellDataSet(as.matrix(count_matrix),#as(count_matrix, "sparseMatrix"),
                                 phenoData = pd,
                                 featureData = fd,
                                 expressionFamily=tobit(),
                                 lowerDetectionLimit=0.1)
rpc_matrix <- relative2abs(EPI.trajectory)
EPI.trajectory <- newCellDataSet(as(as.matrix(rpc_matrix),"sparseMatrix"), phenoData = pd, featureData = fd, expressionFamily=negbinomial.size(), lowerDetectionLimit = 0.5)
EPI.trajectory <- estimateSizeFactors(EPI.trajectory)
EPI.trajectory <- estimateDispersions(EPI.trajectory)

## Genes expressed in at least 1% cells with expressoin > 0.1
EPI.trajectory <- detectGenes(EPI.trajectory, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(EPI.trajectory), num_cells_expressed > nrow(sample_sheet) * 0.005))

disp_table <- dispersionTable(EPI.trajectory)
ordering_genes <- as.character(subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id)

EPI.trajectory <- setOrderingFilter(EPI.trajectory,ordering_genes = ordering_genes) #, 


EPI.trajectory <- reduceDimension(EPI.trajectory,max_components = 2, method = 'DDRTree')
EPI.trajectory <- orderCells(EPI.trajectory,reverse=T)

p1 = plot_cell_trajectory(EPI.trajectory,color_by = "Pseudotime",show_branch_points = FALSE)
p2 = plot_cell_trajectory(EPI.trajectory, color_by = "seurat_clusters",show_branch_points = FALSE,cell_size = 2)+scale_colour_manual(values =color)
pdf("EPI_resubcluster_trajectory.pdf",width = 10,height = 4)
p1+p2
dev.off()
