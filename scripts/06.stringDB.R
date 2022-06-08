library(STRINGdb)
library(org.Hs.eg.db)
library(dplyr)
library(clusterProfiler)
#library(tidyverse)

combine.DE = readRDS("combine.aneuploid-vs-E.GLMM.DEG.result.rds")
combine.DE = combine.DE[which(rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.1 & rowMax(as.matrix(combine.DE[,c("log.mean1","log.mean2")]))>0.1),]
chrom.state = c("E",paste0(c(1:22,"X","Y"),"T"),paste0(c(1:22,"X"),"M"))

combine.DE.filter = combine.DE[which(abs(combine.DE$logFC)>0.25 & combine.DE$p.adj.BH<0.05 & rowMax(as.matrix(combine.DE[,c("pct1","pct2")]))>0.25),]

# creat STRINGdb object
string_db <- STRINGdb$new( version="11.5", species=9606, 
                           score_threshold=400, input_directory="")


string.net = c()
for(i in chrom.state[-1]){
gene = combine.DE.filter[which(combine.DE.filter$chrom.state %in% i),"gene.name"]
gene = data.frame(gene.name = unique(gene),stringsAsFactors=F)

data_mapped <- string_db$map(gene,my_data_frame_id_col_names = "gene.name", 
                                      removeUnmappedRows = TRUE)
info <- string_db$get_interactions(data_mapped$STRING_id)
data_mapped = data_mapped[!duplicated(data_mapped$STRING_id),]
rownames(data_mapped) = data_mapped$STRING_id
info$form.gene = data_mapped[as.character(info$from),"gene.name"]
info$to.gene = data_mapped[as.character(info$to),"gene.name"]
info = info[!duplicated(paste0(info$to,info$from)),]
info$chrom.state = i
string.net = rbind(string.net,info)
}
string.net = string.net[!duplicated(paste0(string.net$form.gene,string.net$to.gene,string.net$chrom.state)),]

write.table(string.net,"each.chrom.state.DEGs.string.network.infor.score400.txt",sep="\t",quote = F,row.names=F)
