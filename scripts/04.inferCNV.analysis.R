library(infercnv)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)
color = c(brewer.pal(9, "Set1"),brewer.pal(8,"Set2")[1:8],brewer.pal(12,"Paired")[1:12],brewer.pal(8,"Dark2")[1:8],brewer.pal(8,"Accent"))
celltype.color = brewer.pal(8,"Set2")[c(4,2,3,1)]
celltype.color = c("#DB6260","#B4B13E","#336AAD","#A13C84")
names(celltype.color) = c("EPI","PE","Polar.TE","Mural.TE")

annot = read.table("meta.data.annot.txt",header=F,sep="\t",row.names=1)
gene.order = read.table("gene.order.txt",header=F,sep="\t",row.names=1)
name = paste0("ENSG",unlist(lapply(rownames(gene.order), function(x) strsplit(x,"-ENSG",fixed = TRUE)[[1]][2])))
gene.name = as.matrix(rownames(gene.order))
rownames(gene.name) = name

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                    annotations_file=annot,
                                    gene_order_file=gene.order,
				    chr_exclude = c("chrM"),
                                    ref_group_names=c("E334","E335","E390","E391","E392","E394"))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="result.infercnv",  # dir is auto-created for storing outputs
                             cluster_by_groups=T,   # cluster
                             denoise=T,
                             HMM=T,
			                       num_threads = 16,
                             window_length = 101,
                             scale_data = T,
)

################heatmap for the result of infercnv#################
infercnv_obj = readRDS("run.final.infercnv_obj")
expr = infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
abnormal_loc <- infercnv_obj@observation_grouped_cell_indices
gn <- rownames(expr)
geneFile <- read.table('E:/project/embryo/data/database/gene.order.txt',stringsAsFactors = F)


color = rep(c("azure4","azure2"),11)
new_cluster = factor(gsub("chr","",geneFile$V2),levels=c(paste0("",c(1:22))))
top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color,col = "white"), # ÉèÖÃÌî³äÉ«
                       labels = levels(new_cluster), #NULL,#
                       labels_gp = gpar(cex = 0.9, col = "white")), show_annotation_name = F) 


ab.type = as.character(gsub("[0-9]","",meta$chrom.statePGS))
ab.type = factor(ab.type,levels = c("E","T","M"))
MTN.color = brewer.pal(8,"Dark2")[1:3]
names(MTN.color) = c("E","T","M")
ha = rowAnnotation(ab.type = as.character(ab.type),
                   celltype = as.character(meta$Lineage),
                   col = list(ab.type = MTN.color,celltype = celltype.color),
                   show_legend = T,
                   #width = unit(1, "cm"),
                   gap = unit(1, "points"),
                   show_annotation_name = F)

mat_scaled = scale(as.matrix(expr),center = T)
mat_scaled[mat_scaled>10] = 10
mat_scaled[mat_scaled< -10] = -10
mat_scaled = mat_scaled/10

col1 = colorRampPalette(c("#b2182b","#f4a582"))(7)
col2 = colorRampPalette(c("#fddbc7","#d1e5f0"))(4)
col3 = colorRampPalette(c("#92c5de", "#2166ac"))(7)
colours = rev(c(col1,col2,col3))

ht_ab = Heatmap(as.matrix(t(mat_scaled)),
                cluster_rows = F,
                cluster_columns = F,
                show_column_names = F,
                show_row_names = F,
              #  col = colours,
                left_annotation = ha,
                top_annotation = top_color,
                show_heatmap_legend=T,
                column_split = new_cluster,
                row_title = "",
                row_title_side = c("right"),
                row_title_rot = 90,
                row_title_gp = gpar(fontsize = 25),
                column_names_gp = gpar(fontsize = 8),
                column_title = NULL, 
                heatmap_legend_param = list(
                  title = " ",
                  title_position = "topcenter",#"leftcenter-rot", 
                  title_gp = gpar(fontsize = 20),
                  at=c(min(mat_scaled),max(mat_scaled)), 
                  legend_height = unit(6, "cm")),
                width = 20, height = 5
) 

pdf("14082cells.NN.ab.sample_chrom.state.rmsex.add.bar.defined.E.scaleTo1.test.pdf",width = 20,height = 25)
draw(ht_ab)
dev.off()

################calculate the difference signal between each chromosome for each cell#########
colnames(geneFile) = c("gene.id","chrom","start","end")
rownames(geneFile) = geneFile$gene.id
infer.obj <- CreateSeuratObject(counts = t(expr), project = "inferCNV.DE", min.cells = 2,meta.data = geneFile, min.features = 100)
infer.obj@meta.data$chrom = as.character(infer.obj@meta.data$chrom)
Idents(infer.obj) = "chrom"
DE = FindAllMarkers(infer.obj, only.pos = F, min.pct = 0.1, logfc.threshold = 0.01,return.thresh = 1)
DE$gene = gsub("-","_",DE$gene)
chrom.DE = merge(DE,meta.data,by.y = "Cell_ID",by.x = "gene")
write.table(chrom.DE,"00.each.chrom.DE.cells.fc0.01.txt",sep = "\t",row.names = F,quote = F)

################summary the result of infercnv#################
chrom.DE = chrom.DE[-which(chrom.DE$cluster %in% c("chrX","chrY")),]
chrom.DE = chrom.DE[-which((chrom.DE$chrom) %in% c("X","Y")),]

chrom.DE.filter = chrom.DE[which(chrom.DE$p_val_adj<0.01 & abs(chrom.DE$avg_logFC)>0.1),]

chrom.DE.filter$sample.infor = as.character(chrom.DE.filter$sample.infor)

sample = unique(meta.data$sample.infor)
sample = as.character(sample[-grep("X|Y",sample)])
ab.number = data.frame()
for(j in 1:length(sample)){
  ab.number[sample[j],"total"] = dim(meta.data[grep(sample[j],meta.data$sample.infor),])[1]
  temp = chrom.DE.filter[which(chrom.DE.filter$sample.infor %in% sample[j]),]
  for(i in 1:22){
    t1 = temp[which((temp$cluster) %in% paste0("chr",i) & temp$avg_logFC>0),]
    ab.number[sample[j],paste0(i,"T")] = dim(t1)[1]
    t2 = temp[which((temp$cluster) %in% paste0("chr",i) & temp$avg_logFC<0),]
    ab.number[sample[j],paste0(i,"M")] = dim(t2)[1]
  }
  ab.number[sample[j],"normal"] = ab.number[sample[j],"total"]-length(unique(as.character(temp$gene)))
}
temp = meta.data[!duplicated(meta.data$sample.infor),]
temp = temp[-which(temp$chrom %in% c("X","Y")),]
temp$chrom = factor(temp$chrom,levels = c(1:22,"E"))
temp$ab.type = factor(temp$ab.type,levels = c("T","M","E"))
temp = temp[order(temp[,"chrom"],temp[,"ab.type"]),]
temp$sample.infor = factor(temp$sample.infor,levels = temp$sample.infor)
ab.number = ab.number[as.matrix(temp$sample.infor),]

write.table(ab.number,"../each_embryo_abtype_cellNumber_in.chrom.txt",sep = "\t",quote = F)