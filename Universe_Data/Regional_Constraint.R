#havrilla CCR
#ccr <- as.data.frame(read.table("Gene_lists/Regional_Constraint/ccrs.autosomes.v2.20180420.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE))
#ccr95 <- ccr[which(ccr$V4>95),]
#ccr95$V5<-checkGeneSymbols(ccr95$V5,unmapped.as.na=FALSE,map=hgnc.table)[[3]]
#save(ccr95, file="output/Data/ccr95.rda", compress="bzip2")
#ccr99 <- ccr[which(ccr$V4>99),]
#ccr99$V5<-checkGeneSymbols(ccr99$V5,unmapped.as.na=FALSE,map=hgnc.table)[[3]]
#save(ccr99, file="output/Data/ccr99.rda", compress="bzip2")

#mult99<-read.table("Gene_lists/Regional_Constraint/genesw99CCR(supp_table_1).tsv", sep = '\t', header = TRUE,stringsAsFactors=FALSE)
#ccr99_genes <- data.frame(gene=unique(ccr99$V5))
#ccr99_genes$n <- ifelse(ccr99_genes$gene%in%mult99$gene,vlookup(ccr99_genes$gene,mult99,result_column="n",lookup_column = "gene"),1)
#save(ccr99_genes, file="output/Data/ccr99_genes.rda", compress="bzip2")
#rm(mult99)
load("output/Data/ccr99_genes.rda")

#highest_ccr <- data.frame(gene=unique(ccr$V5))
#highest_ccr$ccr <- rep(NA,length(highest_ccr$gene))
#highest_ccr$ccr <- unlist(lapply(highest_ccr$gene,function(x) max(ccr$V4[which(ccr$V5==x)])))
#highest_ccr$gene<-checkGeneSymbols(highest_ccr$gene,unmapped.as.na=FALSE,map=hgnc.table)[[3]]
#save(highest_ccr, file="output/Data/highest_ccr.rda", compress="bzip2")
load("output/Data/highest_ccr.rda")

#samocha regional missense constraint
samocha<-read.xlsx("Gene_lists/Regional_Constraint/148353-3.xlsx",sheet=2)
genes_w_var <- data.frame(gene=unique(samocha$gene))
genes_w_var$lowest_gamma <- unlist(lapply(genes_w_var$gene,function(x) min(samocha$obs_exp[which(samocha$gene==x)])))
genes_w_var$gene<-checkGeneSymbols(genes_w_var$gene,unmapped.as.na=FALSE,map=hgnc.table)[[3]]
rm(samocha)

#coban-akdemir NMD- constraint
nmd_min <-read.xlsx("Gene_lists/Regional_Constraint/1-s2.0-S0002929718302039-mmc5.xlsx",startRow=2)
nmd_min<-nmd_min[,c(1,2)]
nmd_min$gene<-checkGeneSymbols(nmd_min$gene,unmapped.as.na=FALSE,map=hgnc.table)[[3]]



