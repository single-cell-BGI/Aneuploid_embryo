library(Seurat)
library(Rmisc)
load("CNVfiltered.IntegrateData.RData")
get.id = function(y){y = paste0("^",y,"-ENSG");grep(pattern = y,x = rownames(embryo.integrated_CNVdefined@assays$RNA@counts),value = T)}
chrom.state = c("E",paste0(c(1:22,"X","Y"),"T"),paste0(c(1:22,"X"),"M"))

term.names = c("Apoptotic","Stress response","Mitochondrial ATP synthesis","Protein modification and metabolism")
names(term.names) = c("M.up.apoptotic","M.up.stress.response","MT.up.mitochondrial","MT.down.protein.folding")
for(term in names(term.names)){
  term.gene = read.table(paste0("select.gene.score/",term,".txt"),sep = "\t",stringsAsFactors = F)
  term.gene = as.list(data.frame(term.gene = unlist(lapply(as.list(term.gene$V1),get.id)),stringsAsFactors = F))
    
  embryo.integrated_CNVdefined <- AddModuleScore(
    object = embryo.integrated_CNVdefined,
    features = term.gene,
    assay = "RNA",
    ctrl = length(unlist(term.gene)), 
    name = paste0(term.names[term],'_score')
  )
    }
