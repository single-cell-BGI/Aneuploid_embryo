########
#library(SCnorm)
library(lmerTest)
library(SingleCellExperiment)
library(data.table)
library(emmeans)
library(dplyr)
library(reshape2)

args = commandArgs(trailingOnly = TRUE)

path = args[1]
norm.sce = args[2]
meta.data = args[3]
chrom.state = args[4]

#####loading normalized data and cell meta infor
norm.sce = readRDS(norm.sce)
meta = read.table(meta.data,sep = "\t",header = T,stringsAsFactors = F)
rownames(meta) = meta$Cell.Id
#norm.sce = results(norm.sce,type = "NormalizedData")

#####subset compared aneuploid and euploid data
meta = meta[which(meta$Karyotype.By.InferCNV %in% c("E",chrom.state)),]

meta$Karyotype.By.InferCNV = factor(meta$Karyotype.By.InferCNV,levels = c("E",chrom.state))
meta$Lineage = factor(meta$Lineage,levels = c("EPI","PE","Polar.TE","Mural.TE"))
norm.sce = norm.sce[,rownames(meta)]

##### data filtering and summary the genes expression infor
means = aggregate(t(norm.sce),by = list(Karyotype.By.InferCNV = meta$Karyotype.By.InferCNV,Lineage = meta$Lineage),FUN=mean)
means = melt(means,c("Karyotype.By.InferCNV","Lineage"),variable.name = "gene.id",value.name = "log.means")
means$log.means = log(means$log.means+1)
means = dcast(means,gene.id+Lineage~Karyotype.By.InferCNV,value.var = "log.means")
colnames(means)[3:4] = c("log.mean2","log.mean1")

pct = aggregate(t(norm.sce),by = list(Karyotype.By.InferCNV = meta$Karyotype.By.InferCNV,Lineage = meta$Lineage),FUN=function(x){round(length(x[x>0])/length(x),2)})
pct = melt(pct,c("Karyotype.By.InferCNV","Lineage"),variable.name = "gene.id",value.name = "pct")
pct = dcast(pct,gene.id+Lineage~Karyotype.By.InferCNV,value.var = "pct")
colnames(pct)[3:4] = c("pct2","pct1")

celltype.infor = merge(pct,means,by = c("gene.id","Lineage"))
celltype.infor$logFC = celltype.infor$log.mean1-celltype.infor$log.mean2
celltype.infor$gene.id = as.character(celltype.infor$gene.id)
########select genes expression at least 10% cells of euploidy or aneuploidy to calculated p value
celltype.infor = celltype.infor[which(celltype.infor$pct2>=0.1|celltype.infor$pct1>=0.1),]
dim(celltype.infor)
#####defined function to calc de
GLMM.celltype.de.f = function(norm.count,meta.temp,gene.id){
  meta.temp$scnorm.counts = as.numeric(norm.count[gene.id,])
  m = glmer.nb(data = meta.temp,formula = round(scnorm.counts + 1) ~ (1 | Embryo.Id) + Karyotype.By.InferCNV*Lineage,nAGQ = 0)
  anova.result = as.data.table(emmeans(m, revpairwise ~ Karyotype.By.InferCNV | Lineage)$contrasts)
  anova.result$gene.id = gene.id
  anova.result$isSingular = isSingular(m)
  return(anova.result)
}


######calculate each genes de p value
library(parallel)
cl <- makeCluster(8) ##Call multithreading: 8 thread ;parallel::detectCores()-1
clusterEvalQ(cl,c(library(lmerTest),library(emmeans),library(data.table))) ####Load R package to multithread
clusterExport(cl,varlist=c("GLMM.celltype.de.f","norm.sce","meta"),env = environment()) ####Load R environment variables to multithread

Sys.time()
DE.celltype = parLapply(cl,as.matrix(unique(celltype.infor$gene.id)),function(x){return(tryCatch(GLMM.celltype.de.f(norm.sce,meta,x), error = function(e) NULL))}) #skip the genes can not converge with error report: Error in f_refitNB pwrssUpdate did not converge in (maxit) iterations
#DE.celltype = lapply(as.matrix(unique(celltype.infor$gene.id)),function(x){print(x);return(tryCatch(GLMM.celltype.de.f(norm.sce,meta,x), error = function(e) NULL))}) #skip the genes can not converge with error report: Error in f_refitNB pwrssUpdate did not converge in (maxit) iterations 
DE.celltype = do.call(rbind,DE.celltype)
Sys.time()

summary.celltype.DE = merge(celltype.infor,DE.celltype,by = c("gene.id","Lineage"),all.x = T)
write.table(summary.celltype.DE,paste0(path,"/",chrom.state,"-vs-E.GLMM.DEG.celltype.specific.result.txt"),quote=F,row.names = F,sep="\t")


