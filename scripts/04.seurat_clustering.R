library(Seurat)

TPM = readRDS("embryo.14082cells.TPM.new.rds")
metadata = read.table("embryo_14082cells.metadata.upload.txt",row.names=F,header=T,stringsAsFactors=F,sep="\t")
rownames(metadata) = metadata$Cell.Id
metadata = metadata[colnames(TPM),]
embryo = CreateSeuratObject(counts = TPM,meta.data = metadata,min.cells = 3, min.features = 800,names.delim = "_")

embryo.list <- SplitObject(embryo, split.by = "Platform")
for (i in 1:length(embryo.list)) {
  embryo.list[[i]] <- NormalizeData(embryo.list[[i]], verbose = FALSE)
  embryo.list[[i]] <- FindVariableFeatures(embryo.list[[i]], selection.method = "vst", 
                                             nfeatures = 20000, verbose = FALSE)
}
rm(embryo,TPM)
embryo.anchors <- FindIntegrationAnchors(object.list = embryo.list, dims = 1:80,anchor.features = 20000)
embryo.integrated <- IntegrateData(anchorset = embryo.anchors, dims = 1:80)
DefaultAssay(embryo.integrated) <- "integrated"
embryo.integrated <- ScaleData(embryo.integrated, verbose = FALSE)
embryo.integrated <- RunPCA(embryo.integrated, npcs = 10, verbose = FALSE)
embryo.integrated <- FindNeighbors(embryo.integrated, dims = 1:10,k.param = 50)
embryo.integrated <- FindClusters(object = embryo.integrated,resolution = 0.2,n.start = 50)
embryo.integrated <- RunUMAP(embryo.integrated, reduction = "pca", dims = 1:10,min.dist = 0.4,n.neighbors = 40,seed.use = 10,negative.sample.rate = 20,set.op.mix.ratio = 0.7) #
embryo.integrated <- RunTSNE(object = embryo.integrated, dims = 1:10, perplexity=100)

DimPlot(embryo.integrated, reduction = "umap", group.by = "Seurat.Clusters",cols= brewer.pal(8,"Pastel2"))
saveRDS(embryo.integrated,"14082cells.IntegrateData.clustering.rds")

Idents(embryo.integrated) = "Lineage"
DefaultAssay(embryo.integrated) = "RNA"
DE.gene = FindAllMarkers(embryo.integrated,only.pos = T, min.pct = 0.25, logfc.threshold = 0.1,return.thresh = 0.05)
DE.gene.filter = DE.gene[DE.gene$avg_logFC>0.25 & DE.gene$p_val_adj<0.05,]
write.table(DE.gene.filter,"14082cells.meta.data.4celltype.FindAllMarkers.txt",sep = "\t",quote = F)

##############summarize the lineage in each aneuploidy and euploidy################
library(Rmisc)
library(ggsignif)
tab.meta = data.frame(table(meta14082$sample.infor,meta14082$Lineage),stringsAsFactors = F)
tab.meta = dcast(tab.meta,Var1~Var2)
tab.meta$total.cell = rowSums(tab.meta[,2:5])

tab.meta.m = melt(tab.meta[,c(2:6,10)],id.vars = "ab.type")
tgc <- summarySE(tab.meta.m, measurevar="value", groupvars=c("variable","ab.type"))
tgc$ab.type = factor(tgc$ab.type,levels = c("E","T","M"))
colnames(tgc)[c(1,4)] = c("celltype","cellnumber.mean")

write.table(tgc,"14082cells.4celltype.NTM.each.embryo.cellnumber.mean.combine.signif.txt",sep = "\t",quote = F,row.names = F)

tab.meta.m = melt(tab.meta[,c(2:6,9)],id.vars = "chrom.state")
tgc <- summarySE(tab.meta.m, measurevar="value", groupvars=c("variable","chrom.state"))
tgc$ab.type = gsub("[0-9]|X|Y","",tgc$chrom.state)
tgc$ab.type = factor(tgc$ab.type,levels = c("E","T","M"))
# tgc$chrom = paste0("chr",gsub("T|M","",tgc$chrom.state))
# tgc$chrom = gsub("chrE","E",tgc$chrom)
# tgc$chrom = factor(tgc$chrom,levels = unique(tgc$chrom))
colnames(tgc)[c(1,3,4)] = c("celltype","embryo.number","cellnumber.mean")

for(i in c("total.cell","EPI","PE","Polar.TE","Mural.TE")){
  for(j in chrom.state){
    annot_1 = wilcox.test(x=tab.meta[tab.meta$chrom.state == "E",i],y=tab.meta[tab.meta$chrom.state == j,i])
    tgc[which(tgc$celltype %in% i & tgc$chrom.state %in% j),"p.value"] = annot_1$p.value
  }
  
}
tgc[tgc$p.value>0.05,"signif"] = ""
tgc[tgc$p.value<0.05,"signif"] = "*"
tgc[tgc$p.value<0.01,"signif"] = "**"
tgc[tgc$p.value<0.001,"signif"] = "***"

p = list()
tgc$ab.type = factor(tgc$ab.type,levels = c("M","T","E"))
for(i in c("total.cell","EPI","PE","Polar.TE","Mural.TE")){
  temp = tgc[tgc$celltype==i,]
#  temp = temp[order(temp$ab.type,temp$cellnumber.mean,decreasing = T),]
#  temp$chrom.state = factor(temp$chrom.state,levels = temp$chrom.state)
  yposition=max(temp$cellnumber.mean+temp$sd,na.rm = T)
  if(i %in% c("total.cell","Mural.TE")){ ylim=max((temp$cellnumber.mean+2*temp$se),na.rm = T)}else{
  ylim=max((temp$cellnumber.mean+1.3*temp$se),na.rm = T)}
  p[[i]] = ggplot(temp, aes(x=chrom.state, y=cellnumber.mean, fill=ab.type)) + 
    geom_bar(position=position_dodge(preserve = "single"), stat="identity",width = 0.6) + #geom_text(mapping = aes(y = cellnumber.mean+se ,label = round(cellnumber.mean,1)),vjust = -1,position = position_dodge(width = .1),size = 2)+
    geom_errorbar(aes(ymin=cellnumber.mean-se, ymax=cellnumber.mean+se),
                  width=.2, # 
                  position=position_dodge(.9,preserve = "single"))+ylab(i)+xlab(NULL)+#ylim(0,max(max(temp$cellnumber.mean+2*temp$se,na.rm = T)))+
    scale_fill_manual(values = MTN.color)+theme_bw()+theme(axis.text.x=element_text(angle=45, hjust=1),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    geom_text(aes(y=!!ylim ,label = signif),position=position_dodge(.9))
  
}
pdf("14082cells.4celltype.NTM.each.embryo.cellnumber.mean.combine.signif.by.chrom.state.pdf",width = 10,height = 10)
plot_grid(plotlist = p,ncol = 1)
dev.off()


tab.meta[,2:5] = tab.meta[,2:5]/rowSums(tab.meta[,2:5])
tab.meta.m = melt(tab.meta[,c(2:6,10)],id.vars = "ab.type")
tgc <- summarySE(tab.meta.m, measurevar="value", groupvars=c("variable","ab.type"))
tgc$ab.type = factor(tgc$ab.type,levels = c("E","T","M"))
colnames(tgc)[c(1,4)] = c("celltype","cellnumber.mean")
tab.meta.m = melt(tab.meta[,c(2:6,9)],id.vars = "chrom.state")
tgc <- summarySE(tab.meta.m, measurevar="value", groupvars=c("variable","chrom.state"))
tgc$ab.type = gsub("[0-9]|X|Y","",tgc$chrom.state)
tgc$ab.type = factor(tgc$ab.type,levels = c("E","T","M"))
colnames(tgc)[c(1,3,4)] = c("celltype","embryo.number","cellnumber.mean")

for(i in c("total.cell","EPI","PE","Polar.TE","Mural.TE")){
  for(j in chrom.state){
    annot_1 = wilcox.test(x=tab.meta[tab.meta$chrom.state == "E",i],y=tab.meta[tab.meta$chrom.state == j,i])
    tgc[which(tgc$celltype %in% i & tgc$chrom.state %in% j),"p.value"] = annot_1$p.value
  }
  
}
tgc[tgc$p.value>0.05,"signif"] = ""
tgc[tgc$p.value<0.05,"signif"] = "*"
tgc[tgc$p.value<0.01,"signif"] = "**"
tgc[tgc$p.value<0.001,"signif"] = "***"

p = list()
for(i in c("EPI","PE","Polar.TE","Mural.TE")){
annot_1<-wilcox.test(x=tab.meta[tab.meta$ab.type == "E",i],y=tab.meta[tab.meta$ab.type == "T",i])
annot_2<-wilcox.test(x=tab.meta[tab.meta$ab.type == "E",i],y=tab.meta[tab.meta$ab.type == "M",i])
annot_3<-wilcox.test(x=tab.meta[tab.meta$ab.type == "T",i],y=tab.meta[tab.meta$ab.type == "M",i])

temp = tgc[tgc$celltype==i,]
yposition=max(temp$cellnumber.mean+temp$sd)
ylim=max(temp$cellnumber.mean+temp$sd)/4

p[[i]] = ggplot(temp, aes(x=ab.type, y=cellnumber.mean, fill=ab.type)) + 
  geom_bar(position=position_dodge(), stat="identity",width = 0.5) + geom_text(mapping = aes(label = round(cellnumber.mean,0)),vjust = -3)+
  geom_errorbar(aes(ymin=cellnumber.mean-se, ymax=cellnumber.mean+se),
                width=.2, #
                position=position_dodge(.9))+xlab(i)+
  scale_fill_manual(values = MTN.color)+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+#facet_wrap(~celltype,ncol = 4)+
geom_signif(annotations = paste0("p=",c(round(annot_2$p.value,4), round(annot_1$p.value,4), round(annot_3$p.value,4))),y_position = c(yposition+ylim,yposition+ylim/2 ,yposition+ ylim/3), xmin=c(0.95, 0.95,2), xmax=c(3, 2,3) )

}

pdf("14082cells.4celltype.ETM.each.embryo.cellnumber.mean.combine.signif.percentage.pdf",width = 12,height = 3.5)
plot_grid(plotlist = p,ncol = 4)
dev.off()
