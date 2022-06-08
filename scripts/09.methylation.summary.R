library(reshape2)
library(dplyr)

chrom.state = c("E",paste0(c(1:22,"X","Y"),"T"),paste0(c(1:22,"X"),"M"))
combine.DE = readRDS("../../combine.aneuploid-vs-E.GLMM.DEG.result.rds")
combine.DE = combine.DE[which(rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.1 & rowMax(as.matrix(combine.DE[,c("log.mean1","log.mean2")]))>0.1),]
combine.DE$DE.type = "undefined"
combine.DE$DE.type[which(abs(combine.DE$logFC)>0.25 & combine.DE$p.adj.BH<0.05 & rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.25)] = "DEGs"
combine.DE$DE.type[which(abs(combine.DE$logFC)<=0.25 & combine.DE$p.adj.BH>0.05 & combine.DE$pct2>0.25)] = "DCGs"
combine.DE$DE.type[which(combine.DE$DE.type %in% "DEGs" & ((combine.DE$logFC>0 & combine.DE$ab.type=="M")|(combine.DE$logFC<0 & combine.DE$ab.type=="T")))] = "DRGs"
cis.combine.DE = combine.DE[which(combine.DE$location %in% "cis"),]


DMR.all = c()
for(i in c("5T","12T","16T","5M","12M","16M")){
  temp = read.table(paste0("05.DMR2gene.melt/E-vs-",i,".merge.intersect.CGmap.sort.DMR.associated_genes.melt.txt"),header = F,stringsAsFactors = F,fill = T,sep = "\t")
  temp$chrom.state = i
  DMR.all = rbind(DMR.all,temp)
}
colnames(DMR.all) = c("chrom","pos.L","pos.R","T.statistics","p.value","E.methylation","ane.methylation","CpG.number","gene.name","basal.gene.list","extended.gene.list","chrom.state")
DMR.all$location = "trans"
DMR.all$location[which(gsub("chr","",DMR.all$chrom)==gsub("M|T","",DMR.all$chrom.state))] = "cis"
DMR.all$up.down = "down"
DMR.all$up.down[grep("-",as.character(DMR.all$T.statistics))] = "up"
###reverse T.statistics
DMR.all$T.statistics[grep("down",DMR.all$up.down)] = -DMR.all$T.statistics[grep("down",DMR.all$up.down)]
DMR.all$T.statistics[grep("up",DMR.all$up.down)] = -DMR.all$T.statistics[grep("up",DMR.all$up.down)]

DMR.all$chrom.state = factor(DMR.all$chrom.state,levels = intersect(chrom.state,unique(DMR.all$chrom.state)))
DMR.all$up.down = factor(DMR.all$up.down,levels = c("up","down"))

DMR.all.cis = DMR.all[which(DMR.all$location %in% "cis"),]
temp = cis.combine.DE[,c("gene.name","chrom.state","DE.type","exp.level")]
DMR.all.cis = merge(DMR.all.cis,temp,by = c("chrom.state","gene.name"),all.x = T)

DMR.all.cis.sub = DMR.all.cis[which(DMR.all.cis$DE.type %in% c("DEGs","DCGs")),]
DMR.all.cis.sub$ab.type = gsub("[0-9]","",DMR.all.cis.sub$chrom.state)
DMR.all.cis.sub$temp = paste0(DMR.all.cis.sub$chrom.state,DMR.all.cis.sub$gene.name)
DMR.all.cis.sub.uniq = DMR.all.cis.sub %>% group_by(temp) %>% top_n(1,-p.value)
DMR.gene.summary = aggregate(DMR.all.cis.sub.uniq$chrom.state,by = list(DMR.all.cis.sub.uniq$chrom.state,DMR.all.cis.sub.uniq$DE.type),FUN = "length")

DMR.all.cis.sub.sig = DMR.all.cis.sub[DMR.all.cis.sub$p.value<0.05,]
DMR.all.cis.sub.sig$temp = paste0(DMR.all.cis.sub.sig$chrom.state,DMR.all.cis.sub.sig$gene.name)
DMR.all.cis.sub.sig.uniq = DMR.all.cis.sub.sig %>% group_by(temp) %>% top_n(1,-p.value)

DMR.all.cis.sub.sig.uniq.used = DMR.all.cis.sub.sig.uniq[grep("^5|12|16",DMR.all.cis.sub.sig.uniq$chrom.state),]
write.table(DMR.all.cis.sub.sig.uniq.used,"aneuploid.5.12.16.TM.DCGs.DEGs.related.cis.DMR.significant.txt",sep = "\t",quote = F,row.names = F)

############对DCG中的T up M down进行热图展示#########
sub.DMR = DMR.all.cis.sub.sig.uniq.used[which(paste0(DMR.all.cis.sub.sig.uniq.used$DE.type,DMR.all.cis.sub.sig.uniq.used$ab.type,DMR.all.cis.sub.sig.uniq.used$up.down) %in% c("DCGsTup","DCGsMdown")),]
write.table(sub.DMR,"aneuploid.5.12.16.TM.DCGs.DEGs.related.cis.DMR.significant.select.272genes.txt",sep = "\t",quote = F,row.names = F)
sub.DMR = sub.DMR[order(sub.DMR$chrom.state,-sub.DMR$T.statistics),]
sub.DMR.T = sub.DMR[which(sub.DMR$ab.type %in% "T"),]
sub.DMR.T = sub.DMR.T[order(sub.DMR.T$chrom.state,-sub.DMR.T$ane.methylation),]
sub.DMR.M = sub.DMR[which(sub.DMR$ab.type %in% "M"),]
sub.DMR.M = sub.DMR.M[order(sub.DMR.M$chrom.state,-sub.DMR.M$E.methylation),]
sub.DMR = rbind(sub.DMR.T,sub.DMR.M)
mtx = as.matrix(sub.DMR[,c("E.methylation","ane.methylation")])
anno = sub.DMR[,c("gene.name","chrom.state")]
library(ComplexHeatmap)
pdf("aneuploid.5.12.16.TM.DCGs.DEGs.related.cis.DMR.significant.methylation.levels1.pdf",width = 2,height = 4)
Heatmap(mtx,row_split = anno$chrom.state,cluster_rows = F,cluster_columns = F,name = "methylation levels",col = viridis_pal()(8)[4:8])
dev.off()
