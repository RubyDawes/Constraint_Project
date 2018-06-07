nmd <- read.xlsx("Gene_lists/NMD_genes/nmd_genes_2018.xlsx")
nmd$gene <- checkGeneSymbols(nmd$gene,unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
nmd$gene[which(nmd$gene=="Sep-09")]<- "SEPT9"
#removing nmd genes that aren't protein coding
nmd <- nmd[-which(!nmd$gene%in%universe&!grepl(pattern="///",nmd$gene)),]
