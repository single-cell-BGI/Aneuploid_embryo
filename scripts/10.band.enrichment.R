
gene.band = read.table("E:/project/embryo/data/database/gencode.v19.chr_patch_hapl_scaff.annotation_gene_info.add.UCSC.hg19.chromsome.band.txt",header = F,sep = "\t",stringsAsFactors = F,fill = T)
gene.band = gene.band[grep("chr",gene.band$V1),]
colnames(gene.band) = c("chr","star","end","ensembl.id","gene.name","band.chr","band.star","band.end","band.id","gieStain")
gene.band$gene.name = gsub("_","-",gene.band$gene.name)
gene.band$gene.id = paste0(gene.band$gene.name,"-",gene.band$ensembl.id)
gene.band$band.id.chr = gsub("chr","",paste0(gene.band$band.chr,gene.band$band.id))

combine.DE = readRDS("E:/project/embryo/data/RNAseq/01.all.embryo.cell/01.filtered_By_CNV/20.GLMM.DEG/00.summary.aneuploid_vs_euploid/combine.aneuploid-vs-E.GLMM.DEG.result.rds")
combine.DE = combine.DE[which(rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.1 & rowMax(as.matrix(combine.DE[,c("log.mean1","log.mean2")]))>0.1),]
combine.DE$DE.type = "undefined"
combine.DE$DE.type[which(abs(combine.DE$logFC)>0.25 & combine.DE$p.adj.BH<0.05 & rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.25)] = "DEGs"
combine.DE$DE.type[which(abs(combine.DE$logFC)<=0.25 & combine.DE$p.adj.BH>0.05 & combine.DE$pct2>0.25)] = "DCGs"
combine.DE$DE.type[which(combine.DE$DE.type %in% "DEGs" & ((combine.DE$logFC>0 & combine.DE$ab.type=="M")|(combine.DE$logFC<0 & combine.DE$ab.type=="T")))] = "DRGs"
cis.combine.DE = combine.DE[which(combine.DE$location %in% "cis"),]
cis.combine.DE.DEGs = cis.combine.DE[which(cis.combine.DE$DE.type %in% "DEGs"),]
cis.combine.DE.DEGs = merge(cis.combine.DE.DEGs,gene.band,by = "gene.id",all.x = T)
###summary the number of expressed genes in each band
gene.band.exp = gene.band[which(gene.band$gene.id %in% unique(cis.combine.DE$gene.id)),]
summary.gene.dist = data.frame(table(gene.band.exp$chr,gene.band.exp$band.id.chr),stringsAsFactors = F)
summary.gene.dist = summary.gene.dist[summary.gene.dist$Freq>0,]
temp = data.frame(table(gene.band.exp$chr),stringsAsFactors = F)
summary.gene.dist = merge(summary.gene.dist,temp,by = "Var1")
colnames(summary.gene.dist) = c("chr","band","band.gene.number","chrom") ###chrom列为整个染色体表达基因数

###summary the number of cis-PDEGs in each band
temp = data.frame(table(cis.combine.DE.DEGs$chr,cis.combine.DE.DEGs$ab.type),stringsAsFactors = F)
temp = dcast(temp,Var1~Var2)
summary.gene.dist.cisDE = merge(summary.gene.dist,temp,by.x = "chr",by.y = "Var1",all.x = T)
temp = data.frame(table(cis.combine.DE.DEGs$band.id.chr,cis.combine.DE.DEGs$ab.type),stringsAsFactors = F)
temp = dcast(temp,Var1~Var2)
summary.gene.dist.cisDE = merge(summary.gene.dist.cisDE,temp,by.x = "band",by.y = "Var1",all.x = T)
colnames(summary.gene.dist.cisDE)[5:8] = c("M.DEgenes.cis.gn","T.DEgenes.cis.gn","M.DEgenes.band.gn","T.DEgenes.band.gn")

###hypergeometric distribution test for each band
N = summary.gene.dist.cisDE$chrom
M = summary.gene.dist.cisDE$band.gene.number
n = summary.gene.dist.cisDE$M.DEgenes.cis.gn
k = summary.gene.dist.cisDE$M.DEgenes.band.gn
summary.gene.dist.cisDE[is.na(summary.gene.dist.cisDE)] = 0
summary.gene.dist.cisDE$M.cis.DE.p = phyper(k-1,M,N-M,n, lower.tail=FALSE)
n = summary.gene.dist.cisDE$T.DEgenes.cis.gn
k = summary.gene.dist.cisDE$T.DEgenes.band.gn
summary.gene.dist.cisDE$T.cis.DE.p = phyper(k-1,M,N-M,n, lower.tail=FALSE)
