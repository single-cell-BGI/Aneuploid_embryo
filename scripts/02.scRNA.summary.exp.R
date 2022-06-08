library(getopt)
commonds = matrix(c("stat","W",1,"character","output","O",1,"character","batch","B",1,"character"),byrow=TRUE, ncol=4)

opt = getopt(commonds)

list.file = list.files(opt$stat, pattern = ".txt")
list.file = list.file[grep("genecount",list.file)]
filelist <- lapply(as.list(paste0(opt$stat,"/",list.file)), read.table, header = T,stringsAsFactors = F, sep = "\t",row.names=1)
count = lapply(filelist,function(x){x = as.data.frame(x[,1])})
count = do.call(cbind,count)
TPM = lapply(filelist,function(x){x = as.data.frame(x[,2])})
TPM = do.call(cbind,TPM)
colnames(count) = gsub("genecount_|.txt","",list.file)
count = cbind("Gene_ID" = rownames(filelist[[1]]),count)

colnames(TPM) = gsub("genecount_|.txt","",list.file)
TPM = round(TPM,3)

gene.anno = read.table("gencode.v19.chr_patch_hapl_scaff.annotation.txt",sep="\t",stringsAsFactors = F,header = T)
rownames(gene.anno) = gene.anno$gene_id
gene.anno = gene.anno[rownames(filelist[[1]])[1:63568],]
gene.add.name = rownames(filelist[[1]])
gene.add.name[1:67049] = gene.anno$gene.id

TPM = cbind("Gene_ID" = gene.add.name,TPM)
count$Gene_ID = TPM$Gene_ID
write.table(count,paste0(opt$output,"/",opt$batch,"RNA_rawcount.xls"),quote=F,sep="\t",row.names=F)
write.table(TPM,paste0(opt$output,"/",opt$batch,"RNA_TPM.xls"),quote=F,sep="\t",row.names=F)
gzip.name = paste0("gzip ",opt$output,"/",opt$batch,"RNA_TPM.xls")
system(gzip.name)
gzip.name = paste0("gzip ",opt$output,"/",opt$batch,"RNA_rawcount.xls")
system(gzip.name)