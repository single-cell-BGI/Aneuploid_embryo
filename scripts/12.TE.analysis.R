library(Seurat)
library(dplyr)
library(monocle)
library(cowplot)
library(ggplot2)
library(viridis)
library(RColorBrewer)
bg<-theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
xbg<-theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=10))
color = c(brewer.pal(12,"Set3")[1:12],brewer.pal(9, "Set1"),brewer.pal(12,"Paired")[1:12],brewer.pal(8,"Dark2")[1:8],brewer.pal(8,"Accent"))


load("embryo.integrated.11561cells.RData")
Polar_Mular_TE<-subset(embryo.integrated,subset=Lineage %in% c("Polar.TE","Mural.TE"))
DefaultAssay(Polar_Mular_TE)<-"RNA"
Variblegene<- FindVariableFeatures(Polar_Mular_TE, selection.method = "vst", nfeatures = 15000)
count_matrix<-Polar_Mular_TE@assays$integrated@data[which(rownames(Polar_Mular_TE@assays$integrated@data) %in% Variblegene),]
	Polar_Mular_TE_trajectory_sample_sheet<-Polar_Mular_TE@meta.data
	Polar_Mular_TE_trajectory_gene_annotation<-data.frame("gene_short_name"=rownames(count_matrix),"gene_name"=rownames(count_matrix))
	row.names(Polar_Mular_TE_trajectory_gene_annotation)<-Polar_Mular_TE_trajectory_gene_annotation$gene_short_name
	pd <- new("AnnotatedDataFrame", data = Polar_Mular_TE_trajectory_sample_sheet)
	fd <- new("AnnotatedDataFrame", data = Polar_Mular_TE_trajectory_gene_annotation)
	Polar_Mular_TE_trajectory <- newCellDataSet(as(count_matrix, "sparseMatrix"),
		phenoData = pd,
                featureData = fd,
                expressionFamily=negbinomial.size(),
		lowerDetectionLimit=1)

	Polar_Mular_TE_trajectory <- estimateSizeFactors(Polar_Mular_TE_trajectory)
	Polar_Mular_TE_trajectory <- estimateDispersions(Polar_Mular_TE_trajectory)
	ordering_genes<-Polar_Mular_TE_trajectory_gene_annotation$gene_short_name
	Polar_Mular_TE_trajectory <- setOrderingFilter(Polar_Mular_TE_trajectory, ordering_genes)
	Polar_Mular_TE_trajectory <- reduceDimension(Polar_Mular_TE_trajectory, max_components = 2,method = 'DDRTree')
	Polar_Mular_TE_trajectory <- orderCells(Polar_Mular_TE_trajectory)

Polar_Mular_TE_trajectory@phenoData@data$phase = pdata$phase
Polar_Mular_TE_trajectory@phenoData@data$ab.type = gsub("N","E",Polar_Mular_TE_trajectory@phenoData@data$ab.type)
Polar_Mular_TE_trajectory@phenoData@data$ab.type = factor(Polar_Mular_TE_trajectory@phenoData@data$ab.type,levels = c("E","T","M"))
plot_cell_trajectory(Polar_Mular_TE_trajectory,color_by = "phase",show_branch_points = FALSE)+scale_colour_manual(values =c(brewer.pal(8,"Dark2")[c(4:7)]))

save(Polar_Mular_TE_trajectory,file="Polar_Mular_TE_trajectory.RData")

##############defined 4 phase################




#############draw heatmap for TE along the pseudotime##########
diff_test_res.pseudo <- differentialGeneTest(Polar_Mular_TE_trajectory,
                                             fullModelFormulaStr = "~sm.ns(Pseudotime, df=4)",
                                             cores = 4) #reducedModelFormulaStr = "~ platform",

diff_test_res.pseudo.filter = diff_test_res.pseudo[diff_test_res.pseudo$qval<0.00001,]
tt = as.matrix(embryo.integrated_CNVdefined@assays$RNA@counts)[as.matrix(diff_test_res.pseudo$gene_short_name),as.matrix(colnames(Polar_Mular_TE_trajectory))]
pp = tt[apply(tt,1,function(x){length(which(x>1))>50}),]
diff_test_res.pseudo.filter = diff_test_res.pseudo.filter[which(diff_test_res.pseudo.filter$gene_short_name %in% rownames(pp)),]

annotation_col = data.frame(cluster = Polar_Mular_TE_trajectory$ab.type,pseudotime = Polar_Mular_TE_trajectory$Pseudotime,stringsAsFactors = F)
annotation_col = annotation_col[order(annotation_col$pseudotime),]
rownames(annotation_col) = 1:length(annotation_col$pseudotime)
newdata <- round(data.frame(Pseudotime = seq(1,length(annotation_col$pseudotime),length.out = 100)),0)
annotation_col = annotation_col[as.matrix(newdata),]
rownames(annotation_col) = 1:length(annotation_col$pseudotime)

pdf("TE.pseudotime.plot.heatmap.pdf",width = 6,height = 8)
p = plot_pseudotime_heatmap(Polar_Mular_TE_trajectory[as.character(diff_test_res.pseudo.filter$gene_short_name),],
                        num_clusters = 4,hclust_method = "ward.D2",
                        cores = 4,return_heatmap = T,
                        trend_formula = "~sm.ns(Pseudotime, df=4)", #sm.ns(Pseudotime, df=3)#split
                        show_rownames = F,add_annotation_col = annotation_col)
dev.off()

##################summarize the distribuction of TE cells along the trajectory############
pdata = pData(Polar_Mular_TE_trajectory)
phase1 = 8
phase2 = 18
phase3 = 30
pdata$phase = "phase2"
pdata[pdata$Pseudotime<=8,"phase"] = "phase1"
pdata[pdata$Pseudotime>=30,"phase"] = "phase4"
pdata[pdata$Pseudotime>18 & pdata$Pseudotime<30,"phase"] = "phase3"


p1 = ggplot(pdata,aes(x=Pseudotime,fill=ab.type))+geom_density(alpha=.7,trim=T)+
  scale_fill_manual(values = MTN.color)+theme_bw()+guides()+
  xlab(NULL)+ylab(NULL)+labs(title = paste0("TE ","distribution"))+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,plot.title = element_text(hjust = 0.5)
  )+geom_vline(xintercept=c(phase1,phase2,phase3), colour="#BB0000", linetype="dashed")#+ylim(0,0.1)

p2 = ggplot(pdata,aes(x=Pseudotime,fill=ChangeNNcelltype))+geom_density(alpha=.7,trim=T)+
  scale_fill_manual(values = celltype.color)+theme_bw()+guides()+
  xlab(NULL)+ylab(NULL)+labs(title = paste0("TE ","distribution"))+
  theme(plot.background = element_blank()
        ,panel.grid.major = element_blank()
        ,panel.grid.minor = element_blank()
        ,panel.border = element_blank()
        ,axis.line = element_line(color = 'black')
        ,plot.title = element_text(hjust = 0.5)
  )+geom_vline(xintercept=c(phase1,phase2,phase3), colour="#BB0000", linetype="dashed")#+ylim(0,0.1)
p1+p2

pdf("TE.pseudotime.MTN.celltype.distribution.phase8.18.30.pdf",width = 14,height = 2)
p1+p2
dev.off()

library(ggbeeswarm)

pdf("04.split4phase/polar.Mural.density.pseudotime.quasirandom.pdf",width = 6,height = 4)
ggplot(as.data.frame(pdata), 
       aes(x = Pseudotime,y = ab.type, 
            colour = ab.type)) +
  geom_quasirandom(groupOnX = F,size = 0.3) +geom_vline(xintercept=c(8,18,30), colour="#BB0000", linetype="dashed")+
  scale_color_manual(values = MTN.color) +  
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + facet_grid(ChangeNNcelltype~.)

dev.off()
