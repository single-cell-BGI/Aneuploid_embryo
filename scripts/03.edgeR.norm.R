raw.counts = readRDS("/da/raid/spwang/project/00.aneuploid.embryo/00.norm.scnorm/embryo.14082cells.raw.counts.rds")
meta14082 = read.table("/da/raid/spwang/project/00.aneuploid.embryo/00.norm.scnorm/embryo_14082cells.metadata.upload.txt",sep = "\t",header = T,stringsAsFactors = F)

rownames(meta14082) = meta14082$Cell.Id
raw.counts = raw.counts[,as.matrix(rownames(meta14082))]

library(edgeR)

filter = apply(raw.counts, 1, function(x) length(x[x>0])>=3)
raw.counts = raw.counts[filter,]
dim(raw.counts)

genes <- rownames(raw.counts)[grep("-ENSG", rownames(raw.counts))]
spikes <- rownames(raw.counts)[grep("^ERCC-", rownames(raw.counts))]

########calc normalized factor by ERCC counts
ERCC.exp = raw.counts[spikes,]
dim(ERCC.exp)
ERCC.y <- DGEList(counts=ERCC.exp, group=meta14082$Karyotype.By.PGT.A)
ERCC.y <- calcNormFactors(ERCC.y)

########
count.norm <- DGEList(counts=raw.counts[genes,], group=meta14082$Karyotype.By.PGT.A)

count.norm$samples$norm.factors = ERCC.y$samples$norm.factors

count.norm <- estimateCommonDisp(count.norm)
#pseudo.counts = ERCC.norm$pseudo.counts

saveRDS(count.norm,"embryo.14082cells.raw.counts.normalized.edgeR.ERCC.rds")
