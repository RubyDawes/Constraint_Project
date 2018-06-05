#get wang, blomen, hart genes and essentialomes
#read wang genes and which are essential in different cell lines
wang <- read.xlsx("Gene_lists/Cell_KOs/wang_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
wang$Gene<- checkGeneSymbols(wang$Gene, unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
wang<- wang[-which(duplicated(wang$Gene)=="TRUE"),] 
wang_essentialome <- unique(wang$Gene[which(wang$essentialome=="essential")])

#read blomen genes and which are essential in different cell lines
blomen <- read.xlsx("Gene_lists/Cell_KOs/blomen_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
blomen$Gene<- checkGeneSymbols(blomen$Gene, unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
blomen<- blomen[-which(duplicated(blomen$Gene)=="TRUE"),] 
blomen_essentialome <- unique(blomen$Gene[which(blomen$essentialome=="essential")])

#read hart genes and which are essential in different cell lines
hart <- read.xlsx("Gene_lists/Cell_KOs/hart_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
hart$Gene<- checkGeneSymbols(hart$Gene, unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
hart<- hart[-which(duplicated(hart$Gene)=="TRUE"),] 
hart_essentialome <- unique(hart$Gene[which(hart$essentialome=="essential")])

common_genes <- intersect(wang$Gene,intersect(hart$Gene,blomen$Gene))
#Seeing essentiality of each gene accross all cell lines in the three papers
cell_KOs <- data.frame(common_genes,ifelse(vlookup(common_genes,blomen,result_column = "hap1",lookup_column="Gene")=="essential",1,0),
                       ifelse(vlookup(common_genes,blomen,result_column = "kbm7",lookup_column="Gene")=="essential",1,0),ifelse(vlookup(common_genes,wang,result_column = "KBM7",lookup_column="Gene")=="essential",1,0),
                       ifelse(vlookup(common_genes,wang,result_column = "k562",lookup_column="Gene")=="essential",1,0),ifelse(vlookup(common_genes,wang,result_column = "jiyoye",lookup_column="Gene")=="essential",1,0),
                       ifelse(vlookup(common_genes,wang,result_column = "raji",lookup_column="Gene")=="essential",1,0),ifelse(vlookup(common_genes,hart,result_column = "hct",lookup_column="Gene")=="essential",1,0),
                       ifelse(vlookup(common_genes,hart,result_column = "hela",lookup_column="Gene")=="essential",1,0),ifelse(vlookup(common_genes,hart,result_column = "gbm",lookup_column="Gene")=="essential",1,0),
                       ifelse(vlookup(common_genes,hart,result_column = "rpe1",lookup_column="Gene")=="essential",1,0),ifelse(vlookup(common_genes,hart,result_column = "dld1",lookup_column="Gene")=="essential",1,0))
colnames(cell_KOs) <- c("Gene","blomen_hap1","blomen_kbm7","wang_kbm7","wang_k562","wang_jiyoye","wang_raji","wang_hct","wang_hela","wang_gbm","wang_rpe1","wang_dld1")
cell_KOs$no_hits <- cell_KOs$blomen_hap1+cell_KOs$blomen_kbm7+cell_KOs$wang_kbm7+cell_KOs$wang_k562+cell_KOs$wang_jiyoye+cell_KOs$wang_raji+
                        cell_KOs$wang_hct+cell_KOs$wang_hela+cell_KOs$wang_hct+cell_KOs$wang_hela+cell_KOs$wang_gbm+cell_KOs$wang_rpe1+cell_KOs$wang_dld1


rm(blomen,hart,wang)
