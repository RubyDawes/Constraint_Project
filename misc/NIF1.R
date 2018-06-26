
load("output/Data/universe_df.rda")
load("output/Data/hgnc.table.rda")

nif1 <- read.xlsx("Gene_lists/NIF1/ncbi_introns_nif1.xlsx")
nif1$gene_id <- checkGeneSymbols(nif1$gene_id,unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]

universe_df$nif1 <- ifelse(universe_df$gene%in%nif1$gene_id,"Y",NA)

save(universe_df, file="output/Data/universe_df.rda", compress="bzip2")
write.xlsx(universe_df,"output/spreadsheets/universe_all_info.xlsx")

length(which(universe_df$omim=="Y"&universe_df$nif1=="Y"))
#195 OMIM genes have bespoke donors

length(which(is.na(universe_df$omim)&universe_df$nif1=="Y"))
#682 non OMIM genes have bespoke donors

