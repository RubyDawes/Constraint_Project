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

#regional constraint
universe_df$highest_ccr<-vlookup(universe_df$gene,highest_ccr,lookup_column = "gene",result_column = "ccr")
universe_df$ccr99<-ifelse(universe_df$gene%in%ccr99_genes$gene,"Y","N")
universe_df$ccr99_n<-vlookup(universe_df$gene,ccr99_genes,lookup_column = "gene",result_column = "n")
rm(ccr99_genes,highest_ccr)

universe_df$regional_missense_constraint <- ifelse(universe_df$gene%in%genes_w_var$gene,"Y","N")
universe_df$regional_missense_gamma<-vlookup(universe_df$gene,genes_w_var,lookup_column = "gene",result_column = "lowest_gamma")
rm(genes_w_var)

universe_df$nmd_min <-ifelse(universe_df$gene%in%nmd_min$gene,"Y","N")
universe_df$nmd_min_rank <- vlookup(universe_df$gene,nmd_min,lookup_column = "gene",result_column = "min.rank")
universe_df$nmd_min_rank<-universe_df$nmd_min_rank/0.01
rm(nmd_min)

universe_df$any_constraint <- ifelse(universe_df$mis_z>=3.09|universe_df$pLI>=0.9|universe_df$ccr99=="Y"|universe_df$regional_missense_constraint=="Y"|universe_df$nmd_min=="Y","Y","N")
universe_df$any_reg_constraint <- ifelse(universe_df$ccr99=="Y"|universe_df$regional_missense_constraint=="Y"|universe_df$nmd_min=="Y","Y","N")

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
universe_df$lethal_IMPC <- vlookup(universe_df$gene,impc,result_column="is_lethal",lookup_column="human_symbol")
universe_df$IMPC_all_MP_ID <- vlookup(universe_df$gene,impc,result_column="all_MP_ID",lookup_column="human_symbol")
universe_df$IMPC_all_MP_phen <- vlookup(universe_df$gene,impc,result_column="all_MP_phen",lookup_column="human_symbol")
universe_df$IMPC_lethal_MP_ID <- vlookup(universe_df$gene,impc,result_column="lethal_MP_ID",lookup_column="human_symbol")
universe_df$IMPC_lethal_MP_phen <- vlookup(universe_df$gene,impc,result_column="lethal_MP_phen",lookup_column="human_symbol")

universe_df$IMPC_ko <- rep(NA,length(universe_df$gene))
universe_df$IMPC_ko[which(lengths(universe_df$IMPC_all_MP_ID)>0)] <-  "Y"

universe_df$mouse_ko <- ifelse(!is.na(universe_df$MGI_ID),"Y",ifelse(universe_df$IMPC_ko=="Y","Y",NA))
universe_df$lethal_mouse <- rep(NA,length(universe_df$gene))
universe_df$lethal_mouse[which(universe_df$mouse_ko=="Y")] <- "N"
universe_df$lethal_mouse[which(universe_df$lethal_MGI=="Y"|universe_df$lethal_IMPC=="Y")] <- "Y"

rm(mgi,impc)

#cell knockouts
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

#human lethal genes- are the genes in my OMIM API search for genes causing lethality in humans
load("output/Data/human_lethal_genes.rda")
universe_df$human_lethal_B <- ifelse(is.na(vlookup(universe_df$gene,lethal_genes[which(lethal_genes$listB=="yes"),],lookup_column = "gene")),"N","Y")
universe_df$human_lethal_A <- ifelse(is.na(vlookup(universe_df$gene,lethal_genes[which(lethal_genes$listA=="yes"),],lookup_column = "gene")),"N","Y")
universe_df$lethal_phen <- lapply(universe_df$gene,function(x) vlookup(x,lethal_genes[which(lethal_genes$listB=="yes"),],lookup_column = "gene",result_column = "lethal_phenotype_mim"))
universe_df$lethal_inheritance <- lapply(universe_df$gene,function(x) vlookup(x,lethal_genes[which(lethal_genes$listB=="yes"),],lookup_column = "gene",result_column = "lethal_inheritance_pattern"))
rm(lethal_genes)

save(universe_df, file="output/Data/universe_df.rda", compress="bzip2")
write.xlsx(universe_df[,c(1:45,48:51)],"output/spreadsheets/universe_all_info.xlsx",append=TRUE)


