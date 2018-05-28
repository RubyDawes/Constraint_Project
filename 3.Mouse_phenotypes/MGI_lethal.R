#plotting proportion of disease, nondisease genes with lethal phenotypes in MGI

source("2.ExAC_constraint/exac_constraint.R")
mgi <- read.xlsx("output/spreadsheets/MGI_genes_with_phenotypes.xlsx")


exac$MGI_ID <- vlookup(exac$gene, mgi,result_column = "MGI_ID",lookup_column = "human_symbol")
exac$MP_terms <- vlookup(exac$gene, mgi,result_column = "MP_ID",lookup_column = "human_symbol")
exac$is_lethal <- vlookup(exac$gene, mgi,result_column = "is_lethal",lookup_column = "human_symbol")


#how many LoF constrained genes?
filter1 <- exac[which(exac$pLI>=0.9),]
#how many of these aren't known disease genes?
filter2 <- filter1[which(filter1$omim=="N"),]
#how many of these have mouse phenotype?
filter3 <- filter2[which(!is.na(filter2$MGI_ID)),]
#how many of these are lethal?
filter4 <- filter3[which(filter3$is_lethal=="Y"),]

  
#how many missense constrained genes?
filter5 <- exac[which(exac$mis_z>=3.09),]
#how many of these aren't known disease genes?
filter6 <- filter5[which(filter5$omim=="N"),]
#how many of these have mouse phenotype?
filter7 <- filter6[which(!is.na(filter6$MGI_ID)),]
#how many of these are lethal?
filter8 <- filter7[which(filter7$is_lethal=="Y"),]
