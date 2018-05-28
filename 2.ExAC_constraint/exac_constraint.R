# a script to take omim genes vs NMD genes vs non-disease genes and compare constraint
source("setup.R")

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
#get omim df
source("1.OMIM_filtering/OMIM_parse_and_filter.R")

#get NMD genes, check symbols
nmd <- read.xlsx("Gene_lists/NMD_genes/nmd_genes_2018.xlsx")
nmd$gene <- checkGeneSymbols(nmd$gene,unmapped.as.na=FALSE)[[3]]

disease_genes <- append(omim$gene ,setdiff(nmd$gene,omim$gene),after=length(omim$gene))

#info on disease status of genes in exac
exac$omim <- ifelse(exac$gene%in%omim$gene,"Y","N")
exac$omim_phen <- ifelse(exac$gene%in%omim$gene,vlookup(exac$gene,omim,result_column="Phenotypes",lookup_column="gene"),NA)
exac$omim_inheritance <- ifelse(exac$gene%in%omim$gene,vlookup(exac$gene,omim,result_column="Inheritance_pattern",lookup_column="gene"),NA)

exac$nmd <- ifelse(exac$gene%in%nmd$gene,"Y","N")
exac$disease <- ifelse(exac$gene%in%disease_genes,"Y","N")



