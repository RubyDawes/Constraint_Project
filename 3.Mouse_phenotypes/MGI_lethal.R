#plotting proportion of disease, nondisease genes with lethal phenotypes in MGI

source("2.ExAC_constraint/exac_constraint.R")
source("3.Mouse_phenotypes/MGI_parse.R")

exac$MGI_ID <- vlookup(exac$gene, mgi,result_column = "MGI_ID",lookup_column = "human_symbol")
exac$MP_terms <- vlookup(exac$gene, mgi,result_column = "MP_ID",lookup_column = "human_symbol")

exac$mouse_lethal <- ifelse(!is.na(vlookup(exac$gene, mgi,result_column = "MGI_ID",lookup_column = "human_symbol")),"Y","N")
