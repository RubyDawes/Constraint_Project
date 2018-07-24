mit <- universe_df[which(grepl(pattern="MT",universe_df$Inheritance_pattern)),]

library("org.Hs.eg.db")
mit$ncbi_id <- unlist(mget(x=as.character(mit$gene),envir=org.Hs.egALIAS2EG))
mit$ncbi_id <-paste("ncbi-geneid:",mit$ncbi_id,sep="")

library(KEGGREST)
pathways<- t(data.frame(as.list(keggList("pathway"))))
kegg_human <- dimnames(t(data.frame(as.list(keggList("hsa")))))[[1]]
kegg_human <-gsub("[.]",":",kegg_human)
kegg <- keggConv("hsa", "ncbi-geneid")
kegg_ncbi <- unlist(lapply(kegg_human,function(x) names(which(kegg==x))))

kegg_id_convert <- data.frame(kegg_human,kegg_ncbi)

mit$kegg_id <- unlist(lapply(mit$ncbi_id, function(x) kegg_id_convert$kegg_human[which(kegg_id_convert$kegg_ncbi==x)]))

mit$CDS <- unlist(lapply(mit$kegg_id,function(x) keggGet(x)[[1]][1]$ENTRY))

mit$kegg_pathways <- lapply(mit$kegg_id, function(x) 
  if ("PATHWAY"%in%names(keggGet(x)[[1]])){keggGet(x)[[1]][which(names(keggGet(x)[[1]])=="PATHWAY")]$PATHWAY} else {"no pathway"})

which(grepl("Oxidative phosphorylation",mit$kegg_pathways))
mit$cell_essential[which(grepl("Oxidative phosphorylation",mit$kegg_pathways))]
checkitout<-mit[which(grepl("Oxidative phosphorylation",mit$kegg_pathways)),c(2,3,4,6,8,9,12,14,28,31,40)]
checkitout$causal_variant <- vlookup(checkitout$gene,universe_df,lookup_column = "gene",result_column = "phenotype_causal_variants_simp")

mit$complex1 <- lapply(1:length(mit$go_names),function(x) ifelse("mitochondrial respiratory chain complex I assembly"%in%unlist(mit$go_names[[x]]),"Complex I",NA))

ar <- universe_df[which(universe_df$Inheritance_pattern=="AR"),]
write.xlsx(ar[1:31],"ar.xlsx")

