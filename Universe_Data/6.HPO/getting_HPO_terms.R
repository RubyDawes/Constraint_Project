
HPO<-read.table("Gene_lists/HPO/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt",header=FALSE,sep="\t",quote="",fill=TRUE,row.names=NULL)
colnames(HPO) <- c("gene_id","gene","HPO_name","HPO_id")


#update gene symbols
HPO$gene <- checkGeneSymbols(HPO$gene,unmapped.as.na = FALSE)[[3]]

universe_df$HPO_id <- lapply(universe_df$gene,FUN=function(x) as.character(unique(HPO$HPO_id[which(HPO$gene==x)])))
universe_df$HPO_id[which(lengths(universe_df$HPO_id)==0)] <- list(NULL)
universe_df$HPO_name <- lapply(universe_df$gene,FUN=function(x) as.character(unique(HPO$HPO_name[which(HPO$gene==x)])))
universe_df$HPO_name[which(lengths(universe_df$HPO_name)==0)] <- list(NULL)

#hpo <- get_ontology("Gene_lists/HPO/hp.obo", propagate_relationships = "is_a",extract_tags = "minimal")

rm(HPO)

hpo_annotations <- data.frame(universe_df$gene,I(universe_df$HPO_id),I(universe_df$HPO_name))
save(hpo_annotations, file="output/Data/HPO_annotations.rda", compress="bzip2")
rm(hpo_annotations)