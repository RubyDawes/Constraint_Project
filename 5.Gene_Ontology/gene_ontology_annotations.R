
goa <- read.table("Gene_lists/Gene_Ontology/goa_human.gaf",header=FALSE,sep="\t",comment.char="!",fill=TRUE)
goa$V7 <- checkGeneSymbols(goa$V7,hgnc.table=hgnc.table,unmapped.as.na = FALSE)[[3]]
b<- lapply(universe_df$gene, function(x) unique(goa$V5[which(goa$V3==as.character(x))]))


write.xlsx(universe_df,"output/spreadsheets/universe_all_info.xlsx")