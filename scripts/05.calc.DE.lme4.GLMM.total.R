########
#library(SCnorm)
library(lmerTest)
library(SingleCellExperiment)
library(data.table)
library(emmeans)
library(dplyr)

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
norm.sce = norm.sce[,rownames(meta)]

##### data filtering and summary the genes expression infor
E.gene = norm.sce[,rownames(meta[which(meta$Karyotype.By.InferCNV %in% "E"),])]
E.gene.infor = data.frame(pct2 = apply(E.gene,1,function(x){round(length(x[x>0])/length(x),2)}),log.mean2 = log(rowMeans(E.gene)+1)) ###calculated pct and log.means of genes expression in euploidy
E.gene = apply(E.gene, 1, function(x) length(x[x>0])>=3) 
ane.gene = norm.sce[,rownames(meta[which(meta$Karyotype.By.InferCNV %in% chrom.state),])]
ane.gene.infor = data.frame(pct1 = apply(ane.gene,1,function(x){round(length(x[x>0])/length(x),2)}),log.mean1 = log(rowMeans(ane.gene)+1)) ###calculated pct and log.means of genes expression in aneuploidy
ane.gene = apply(ane.gene, 1, function(x) length(x[x>0])>=3)
gene.used = union(names(E.gene[which(E.gene)]),names(ane.gene[which(ane.gene)])) #####select genes expression at least 3 cells of euploidy or aneuploidy

bind.infor = cbind(ane.gene.infor,E.gene.infor)
bind.infor$logFC = bind.infor$log.mean1-bind.infor$log.mean2
bind.infor = bind.infor[gene.used,]

########select genes expression at least 10% cells of euploidy or aneuploidy to calculated p value
sub.bind.infor = bind.infor[which(rowMax(as.matrix(bind.infor[,c("pct1","pct2")]))>=0.1),]
dim(sub.bind.infor)
#####defined function to calc de
  GLMM.de.f = function(norm.count,meta.temp,gene.id){
    meta.temp$scnorm.counts = as.numeric(norm.count[gene.id,])
    m1 =  suppressWarnings(glmer.nb(data = meta.temp,formula = round(scnorm.counts + 1) ~ (1 | Embryo.Id)+(1 | Lineage)+Karyotype.By.InferCNV,nAGQ = 0))
    m2 =  suppressWarnings(glmer.nb(data = meta.temp,formula = round(scnorm.counts + 1) ~ (1 + Karyotype.By.InferCNV | Embryo.Id)+(1 + Karyotype.By.InferCNV | Lineage)+ Karyotype.By.InferCNV,nAGQ = 0))
    anova.test <- anova(m2, m1, test = "Chisq")
    if(anova.test[2,]$`Pr(>Chisq)`<0.05){
      anova.result = as.data.table(emmeans(m2, revpairwise ~ Karyotype.By.InferCNV)$contrasts)[,c("contrast","estimate","SE","p.value")]
      anova.result = cbind(gene.id = gene.id,model = "m2",isSingular = isSingular(m2),anova.result)
    }else{
      anova.result = as.data.table(emmeans(m1, revpairwise ~ Karyotype.By.InferCNV)$contrasts)[,c("contrast","estimate","SE","p.value")]
      anova.result = cbind(gene.id = gene.id,model = "m1",isSingular = isSingular(m1),anova.result)
    }
    return(as.data.table(anova.result,drop=F))
  }

######calculate each genes de p value
library(parallel)
cl <- makeCluster(parallel::detectCores()-1) ##Call multithreading
clusterEvalQ(cl,c(library(lmerTest),library(emmeans),library(data.table))) ####Load R package to multithread
clusterExport(cl,varlist=c("GLMM.de.f","norm.sce","meta"),env = environment()) ####Load R environment variables to multithread

Sys.time()
emmeans.result = parLapply(cl,as.matrix(rownames(sub.bind.infor)),function(x){print(x);return(tryCatch(GLMM.de.f(norm.sce,meta,x), error = function(e) NULL))})
emmeans.result = do.call(rbind,emmeans.result)
Sys.time()

de.summary = merge(bind.infor,emmeans.result,by.x = "row.names",by.y = "gene.id",all.x = T)
write.table(de.summary,paste0(path,"/",chrom.state,"-vs-E.GLMM.DEG.result.txt"),quote=F,row.names = F,sep="\t")


