
#source("2.ExAC_constraint/exac_constraint.R") #(calls setup.R,OMIM_parse_and_filter.R)
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
rm(exac)

#MGI mouse knockouts- lethal or non-lethal in a mouse
universe_df$mouse_ko <- ifelse(universe_df$gene%in%mgi$human_symbol,"Y",NA)
universe_df$MGI_ID <- vlookup(universe_df$gene,mgi,result_column="MGI_ID",lookup_column="human_symbol")
universe_df$lethal_mouse <- vlookup(universe_df$gene,mgi,result_column="is_lethal",lookup_column="human_symbol")
universe_df$MP_ID <- vlookup(universe_df$gene,mgi,result_column="MP_ID",lookup_column="human_symbol")
universe_df$MP_phen <- vlookup(universe_df$gene,mgi,result_column="MP_phen",lookup_column="human_symbol")
universe_df$allele_info <- vlookup(universe_df$gene,mgi,result_column="allele_info",lookup_column="human_symbol")
universe_df$mouse_symbol <- vlookup(universe_df$gene,mgi,result_column="mouse_symbol",lookup_column="human_symbol")
rm(mgi)

#cell knockouts
#source("4.Cell_Knockouts/cell_essential_genes.R")
universe_df$cell_ko <- ifelse(universe_df$gene%in%cell_KOs$Gene,"Y",NA)
universe_df$cell_essential_hits <- vlookup(universe_df$gene,cell_KOs,result_column="no_hits",lookup_column="Gene")
rm(cell_KOs)

save(universe_df, file="output/Data/universe_df.rda", compress="bzip2")
write.xlsx(universe_df,"output/spreadsheets/universe_all_info.xlsx")

