#how many genes in common between the three papers' essentialomes?
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

common_essentialome <- intersect(hart_essentialome,intersect(blomen_essentialome,wang_essentialome))
semi_essentialome <- wang_essentialome[which(wang_essentialome%in%blomen_essentialome|wang_essentialome%in%hart_essentialome)]
one_essentialome <- union(hart_essentialome,union(blomen_essentialome,wang_essentialome))

