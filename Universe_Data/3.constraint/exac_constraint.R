# a script to get exac scores for protein coding genes
#getting exac data - now updated to gnomad constraint scores (release 2.1, downloaded 22.1.19)
#taking only needed columns
exac <- read.table("Gene_lists/ExAC/constraint.txt.bgz",header = TRUE,fill = TRUE)


#keeping only scores calculated on canonical transcripts with no gene issues
exac <- exac[which(exac$canonical=="true"&exac$gene_issues=="[]"),c(1,2,20,21,22)]

#removing exac genes that aren't protein coding
exac <- exac[-which(!exac$gene%in%universe_df$gene),]

#remove duplicates
exac <- exac[-which(duplicated(exac$gene)==TRUE),]


# keeping constraint scores from release 2.0 just in case
exac_old <- read.table("Gene_lists/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt",header = TRUE,fill = TRUE)

exac_old <- exac_old[,c(2,17,18,20)]
#checking gene names
genecheck <- checkGeneSymbols(exac_old$gene,unmapped.as.na=FALSE)
baddups <- genecheck[which(duplicated(genecheck[[3]])==TRUE),3]
a<- genecheck[which(genecheck[[3]]%in%baddups),]
b<- which(genecheck[[1]]%in%a[which(a[[2]]==FALSE),1])
genecheck <- genecheck[-b,]
exac_old <-exac_old[-b,]
exac_old$gene <- genecheck[[3]]
rm(a,b,baddups,genecheck)

#removing exac genes that aren't protein coding
exac_old <- exac_old[-which(!exac_old$gene%in%universe_df$gene&!grepl(pattern="///",exac_old$gene)),]



