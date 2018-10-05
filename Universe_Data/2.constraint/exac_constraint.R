# a script to get exac scores for protein coding genes
#getting exac data
#taking only needed columns
exac <- read.table("Gene_lists/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt",header = TRUE,fill = TRUE)
exac <- exac[,c(2,17,18,20)]
#checking gene names
genecheck <- checkGeneSymbols(exac$gene,unmapped.as.na=FALSE)
baddups <- genecheck[which(duplicated(genecheck[[3]])==TRUE),3]
a<- genecheck[which(genecheck[[3]]%in%baddups),]
b<- which(genecheck[[1]]%in%a[which(a[[2]]==FALSE),1])
genecheck <- genecheck[-b,]
exac <-exac[-b,]
exac$gene <- genecheck[[3]]
rm(a,b,baddups,genecheck)

#removing exac genes that aren't protein coding
exac <- exac[-which(!exac$gene%in%universe&!grepl(pattern="///",exac$gene)),]




