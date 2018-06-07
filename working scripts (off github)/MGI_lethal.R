#plotting proportion of disease, nondisease genes with lethal phenotypes in MGI

source("2.ExAC_constraint/exac_constraint.R")
source("plot_functions.R")
mgi <- read.xlsx("output/spreadsheets/MGI_genes_with_phenotypes.xlsx")


exac$MGI_ID <- vlookup(exac$gene, mgi,result_column = "MGI_ID",lookup_column = "human_symbol")
exac$MP_terms <- vlookup(exac$gene, mgi,result_column = "MP_ID",lookup_column = "human_symbol")
exac$is_lethal <- vlookup(exac$gene, mgi,result_column = "is_lethal",lookup_column = "human_symbol")


exac_with_phen <- exac[-which(is.na(exac$MGI_ID)),]
exac_nondisease <- exac_with_phen[-which(exac_with_phen$omim=="Y"),]
exac_disease <- exac_with_phen[which(exac_with_phen$omim=="Y"),]
#how many genes with phens are lethal vs non-lethal?


#how many non-disease genes are lethal vs non-lethal?
omim_nondisease_lethal <- data.frame(group=c("Lethal","Non-lethal"),value=c(length(which(exac_nondisease$is_lethal=="Y")),length(which(exac_nondisease$is_lethal=="N"))))
a <- ggplot(omim_nondisease_lethal, aes(x="", y=value, fill=group))+ggtitle(paste("Non-OMIM genes with mouse phenotype \n n=",length(exac_nondisease$gene)))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
a<- a+scale_fill_manual(values=c("brown3","dodgerblue3"))+geom_text(aes(y = value,label = percent(value/sum(value)),family="Avenir"),size=5,position = position_stack(vjust = 0.5))
#how many omim genes are lethal vs non-lethal
omim_disease_lethal <- data.frame(group=c("Lethal","Non-lethal"),value=c(length(which(exac_disease$is_lethal=="Y")),length(which(exac_disease$is_lethal=="N"))))
b <- ggplot(omim_disease_lethal, aes(x="", y=value, fill=group))+ggtitle(paste("OMIM genes with mouse phenotype \n n=",length(exac_disease$gene)))
b<- b+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
b<- b+scale_fill_manual(values=c("brown3","dodgerblue3"))+geom_text(aes(y = value,label = percent(value/sum(value)),family="Avenir"),size=5,position = position_stack(vjust = 0.5))



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


#universe disease genes with lethal mouse phenotype
g<- universe_df$gene[which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"&universe_df$lethal_mouse=="Y")]
h<- exac_disease$gene[which(exac_disease$is_lethal=="Y")]

