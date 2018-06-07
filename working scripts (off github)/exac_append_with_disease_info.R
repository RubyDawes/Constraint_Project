
#get omim df
source("1.OMIM_filtering/OMIM_parse_and_filter.R")

#get NMD genes, check symbols

disease_genes <- append(omim$Gene ,setdiff(nmd$gene,omim$Gene),after=length(omim$Gene))

#info on disease status of genes in exac
exac$omim <- ifelse(exac$gene%in%omim$Gene,"Y","N")
exac$omim_phen <- ifelse(exac$gene%in%omim$Gene,vlookup(exac$gene,omim,result_column="Phenotypes",lookup_column="Gene"),NA)
exac$omim_inheritance <- ifelse(exac$gene%in%omim$Gene,vlookup(exac$gene,omim,result_column="Inheritance_pattern",lookup_column="Gene"),NA)

exac$nmd <- ifelse(exac$gene%in%nmd$gene,"Y","N")
exac$disease <- ifelse(exac$gene%in%disease_genes,"Y","N")
