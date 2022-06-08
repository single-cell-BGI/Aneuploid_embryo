# Usage: Rscript 01.count_GeneReadNumber_ERCC.R gencode.v19.chr_patch_hapl_scaff.annotation.gtf xxx/bamfile/*.bam gencode.v19.chr_patch_hapl_scaff.annotation.gene.length.txt output.genecount.txt

args= commandArgs(T)

library(GenomicFeatures)
library(GenomicAlignments)
library(data.table)

gtf=args[1]
bam=args[2]
geneLength=args[3]
outfile=args[4]
ercc=args[5]

temp = fread(ercc)

if(dim(temp)[1]>0){table_ercc=as.data.frame(table(fread(ercc,header=F)))
rownames(table_ercc)=table_ercc[,1]
table_ercc=table_ercc[2]
}else{table_ercc = data.frame()
table_ercc["ERCC-00002_DQ459430","Freq"]= 0
}

txdb=makeTxDbFromGFF(gtf)
genes=exonsBy(txdb,"gene")
ex=summarizeOverlaps(genes,bam,mode="Union",inter.feature=TRUE,fragment=TRUE,singleEnd=FALSE)
re=as.data.frame(assay(ex))
colnames(table_ercc)=colnames(re)
re_ercc=rbind(table_ercc,re)
re_ercc$GeneID=rownames(re_ercc)

GeneLength = read.table(geneLength,col.names=c("GeneID","GeneLength"),sep="\t")
merge=merge(GeneLength,re_ercc,by="GeneID",all=T)
merge[is.na(merge)]=0
merge.ERCC = merge[grep("ERCC",merge$GeneID),]
merge.ERCC$TPM = merge.ERCC[,3]
merge = merge[-grep("ERCC",merge$GeneID),]
Read_normalize = as.numeric(as.vector(merge[,3])) / as.numeric(as.vector(merge[,2]))
merge$TPM = round(Read_normalize / sum(Read_normalize) *1000000,3)
merge = rbind(merge,merge.ERCC)
rownames(merge)=merge[,1]
merge=merge[,-c(1,2)]

write.table(merge,file=outfile,quote=FALSE,sep="\t")
