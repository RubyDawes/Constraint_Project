rm(universe)
#disease status, phenotype, inheritance pattern
universe_df$omim <- ifelse(universe_df$gene%in%omim$Gene,"Y",NA)
universe_df$mim_number <- vlookup(universe_df$gene,omim,result_column="Mim.Number",lookup_column="Gene")
universe_df$phenotype <- vlookup(universe_df$gene,omim,result_column="Phenotypes",lookup_column="Gene")
universe_df$Inheritance_pattern <- vlookup(universe_df$gene,omim,result_column="Inheritance_pattern",lookup_column="Gene")

universe_df$nmd <- ifelse(universe_df$gene%in%nmd$gene,"Y",NA)

rm(omim,nmd)
#exac scores
universe_df$exac <- ifelse(universe_df$gene%in%exac$gene,"Y",NA)
universe_df$mis_z <- vlookup(universe_df$gene,exac,result_column="mis_z",lookup_column="gene")
universe_df$syn_z <- vlookup(universe_df$gene,exac,result_column="syn_z",lookup_column="gene")
universe_df$pLI <- vlookup(universe_df$gene,exac,result_column="pLI",lookup_column="gene")
universe_df$constrained <- ifelse(universe_df$mis_z>=3.09|universe_df$pLI>=0.9,"Y","N")
rm(exac)

#MGI mouse knockouts- lethal or non-lethal in a mouse
universe_df$MGI_ID <- vlookup(universe_df$gene,mgi,result_column="MGI_ID",lookup_column="human_symbol")
universe_df$lethal_MGI <- vlookup(universe_df$gene,mgi,result_column="is_lethal",lookup_column="human_symbol")
universe_df$lethal_MP_ID <- vlookup(universe_df$gene,mgi,result_column="lethal_MP_ID",lookup_column="human_symbol")
universe_df$lethal_MP_phen <- vlookup(universe_df$gene,mgi,result_column="lethal_MP_phen",lookup_column="human_symbol")
universe_df$all_MP_ID <- vlookup(universe_df$gene,mgi,result_column="all_MP_ID",lookup_column="human_symbol")
universe_df$all_MP_phen <- vlookup(universe_df$gene,mgi,result_column="all_MP_phen",lookup_column="human_symbol")
universe_df$high_MP_ID <- vlookup(universe_df$gene,mgi,result_column="high_MP_ID",lookup_column="human_symbol")
universe_df$high_MP_phen <- vlookup(universe_df$gene,mgi,result_column="high_MP_phen",lookup_column="human_symbol")
universe_df$allele_info <- vlookup(universe_df$gene,mgi,result_column="allele_info",lookup_column="human_symbol")
universe_df$mouse_symbol <- vlookup(universe_df$gene,mgi,result_column="mouse_symbol",lookup_column="human_symbol")

#IMPC phenotype data
universe_df$IMPC_phen <- vlookup(universe_df$gene,impc,result_column="IMPC",lookup_column="Gene_symbol")

universe_df$mouse_ko <- ifelse(!is.na(universe_df$MGI_ID),"Y",ifelse(!is.na(universe_df$IMPC_phen),"Y",NA))
universe_df$lethal_mouse <- ifelse(grepl(pattern="Y|Lethal",paste(universe_df$lethal_MGI,universe_df$IMPC_phen)),"Y",ifelse(!is.na(universe_df$mouse_ko),"N",NA))
rm(mgi,impc)

#cell knockouts
#source("4.Cell_Knockouts/cell_essential_genes.R")
universe_df$cell_ko <- ifelse(universe_df$gene%in%cell_KOs$Gene,"Y",NA)
universe_df$cell_essential_hits <- vlookup(universe_df$gene,cell_KOs,result_column="no_hits",lookup_column="Gene")
universe_df$cell_essential <- ifelse(universe_df$cell_essential_hits>=3,"Y","N")
rm(cell_KOs)


#GO terms
##slow- don't recommend running
#source("5.GO/gene_ontology_annotations.R)
load("output/Data/GO_annotations.rda")
universe_df$go_terms <- go_annotations$universe_df.go_terms
universe_df$go_names <- go_annotations$universe_df.go_names
rm(go_annotations)

#HPO terms
#this is slow- just load
#source("5.HPO/getting_HPO_terms.R")
load("output/Data/HPO_annotations.rda")
universe_df$hpo_terms <- as.character(hpo_annotations$universe_df.HPO_id)
universe_df$hpo_names <- hpo_annotations$universe_df.HPO_name
rm(hpo_annotations)


save(universe_df, file="output/Data/universe_df.rda", compress="bzip2")
write.xlsx(universe_df[,c(1:31)],"output/spreadsheets/universe_all_info.xlsx",append=TRUE)

