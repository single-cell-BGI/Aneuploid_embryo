library(RIdeogram)
library(rCGH)
library(Biobase)

human_karyotype.hg19 = hg19
human_karyotype.hg19$star = 0
human_karyotype.hg19 = human_karyotype.hg19[,c(1,6,2,3,4)]
colnames(human_karyotype.hg19) = c("Chr","Start","End","CE_start","CE_end")
human_karyotype.hg19[23:24,1] = c("X","Y")
human_karyotype.hg19$Chr = paste0("chr",human_karyotype.hg19$Chr)

gene.order =  read.table('E:/project/embryo/data/database/gene.order.txt',stringsAsFactors = F)
combine.DE = readRDS("E:/project/embryo/data/RNAseq/01.all.embryo.cell/01.filtered_By_CNV/20.GLMM.DEG/00.summary.aneuploid_vs_euploid/combine.aneuploid-vs-E.GLMM.DEG.result.rds")
combine.DE = combine.DE[which(rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.1 & rowMax(as.matrix(combine.DE[,c("log.mean1","log.mean2")]))>0.1),]
combine.DE$DE.type = "undefined"
combine.DE$DE.type[which(abs(combine.DE$logFC)>0.25 & combine.DE$p.adj.BH<0.05 & rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.25)] = "DEGs"
combine.DE$DE.type[which(abs(combine.DE$logFC)<=0.25 & combine.DE$p.adj.BH>0.05 & combine.DE$pct2>0.25)] = "DCGs"
combine.DE$DE.type[which(combine.DE$DE.type %in% "DEGs" & ((combine.DE$logFC>0 & combine.DE$ab.type=="M")|(combine.DE$logFC<0 & combine.DE$ab.type=="T")))] = "DRGs"
cis.combine.DE = combine.DE[which(combine.DE$location %in% "cis"),]
cis.combine.DE = cis.combine.DE[-which(cis.combine.DE$DE.type %in% c("DRGs","undefined","DCGs")),]
cis.combine.DE = merge(cis.combine.DE,gene.order,by.x = "gene.id",by.y = "V1",all.x = T)

#######################T-vs-M overlapd dist#########
TF.list =  read.table('E:/project/embryo/data/database/animalTFDB3.0.1639_human.txt',stringsAsFactors = F)
cis.combine.DE = cis.combine.DE[-which(cis.combine.DE$DE.type %in% c("DRGs","undefined","DCGs")),]
cis.combine.DE = merge(cis.combine.DE,gene.order,by.x = "gene.id",by.y = "V1",all.x = T)
cis.combine.DE$TF[which(cis.combine.DE$gene.id %in% TF.list$V1)] = "TF"
cis.combine.DE$temp = cis.combine.DE$ab.type
cis.combine.DE$temp[duplicated(cis.combine.DE$gene.id)] = "overlap"
cis.combine.DE$temp[which(cis.combine.DE$gene.id %in% cis.combine.DE[cis.combine.DE$temp=="overlap","gene.id"])] = "overlap"
cis.combine.DE = cis.combine.DE[!duplicated(paste0(cis.combine.DE$gene.id,cis.combine.DE$temp)),]

TM.density = cis.combine.DE[,c(25:27,29)]
TM.density$temp[TM.density$temp=="T"] = 1
TM.density$temp[TM.density$temp=="M"] = -1
TM.density$temp[TM.density$temp=="overlap"] = 0
TM.density$temp = as.numeric(TM.density$temp)
colnames(TM.density) = c("Chr","Start","End","Value")

TF.label = cis.combine.DE[which(cis.combine.DE$TF %in% "TF"),]
colnames(TF.label)[25:27] = c("chr","chr.star","chr.end")
write.table(TF.label,"cis.TF.DEGs.list.txt",sep = "\t",quote = F,row.names = F)
TF.label$Shape = "triangle"
TF.label$color = "00ffff"
TF.label = TF.label[,c(28,30,25:27,31)]
colnames(TF.label) = c("Type","Shape","Chr","Start","End","color")
band.label = cis.sig.band[,c(6,5,2:4)]
colnames(band.label)[2] = "Shape"
band.label$Shape = "box"
band.label$color = "00ffbb"
colnames(band.label) = c("Type","Shape","Chr","Start","End","color")
ideogram(karyotype = human_karyotype.hg19, overlaid = TM.density, label = TF.label, label_type = "marker", colorset1 = c( "#386CB0","#F0027F","#BF5B17")) #use the arguments 'colorset1' and 'colorset2' to set the colors for gene and LTR heatmaps, separately.
convertSVG("chromosome.svg", device = "png")
file.rename("chromosome.svg",paste0("ideogram.chromosome.cis.PDEGs.","ideogram",".svg"))
file.rename("chromosome.png",paste0("ideogram.chromosome.cis.PDEGs.","ideogram",".png"))



