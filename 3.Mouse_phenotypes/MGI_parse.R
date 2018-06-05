#a script to parse MGI_PhenoGenoMP to extract all genes with lethal phenotypes, write to xlsx
#load in all genes and phenotypes

mgi <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenoGenoMP.xlsx")
mgi <- mgi[,c(1,2,6,4)]
mgi<- mgi[which(!grepl(pattern = "\\,", x = mgi$MGI_ID, ignore.case = FALSE)),]

allele <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenotypicAllele.xlsx")
mgi$allele_info <- vlookup(mgi$`Allele.Symbol(s)`,allele,result_column="Allele.Attribute",lookup_column="Allele.Symbol")
mgi <- mgi[which(grepl(pattern="Null/knockout",x=mgi$allele_info,ignore.case=FALSE)),]

#which genes have a lethal phenotype
mgi_lethalphens <- read.xlsx("Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx")
mgi$is_lethal <- ifelse(mgi$MP_ID%in%mgi_lethalphens$MP.id,"Y","N")
mgi$MP_ID[which(mgi$is_lethal=="N")]<- NA
#concatenating lethal phenotypes so they're all in one row
lethal_genes <- unique(mgi$MGI_ID[which(mgi$is_lethal=="Y")])

for (p in 1:length(lethal_genes)) {
  rows <- which(mgi$MGI_ID==lethal_genes[p])
  a<- unique(mgi$MP_ID[rows])
  mgi$MP_ID[rows[1]] <- paste(a[which(!is.na(a))],collapse=", ")
  if (length(rows)>1) {  
    mgi <- mgi[-rows[-1],]
  }
  print(p)
}

#get rid of non-lethal gene duplicate rows
mgi$is_lethal[which(!is.na(mgi$MP_ID))]<-"Y"
mgi<- mgi[-which(duplicated(mgi$MGI_ID)),]

# add on gene names and human ortholog names
mgi_names <- read.xlsx("Gene_lists/MGI_MPs/HMD_HumanPhenotype.xlsx")
mgi_names$MGI.Marker.Accession.ID <- gsub(" ", "", mgi_names$MGI.Marker.Accession.ID, fixed = TRUE)
mgi$mouse_symbol <- vlookup(mgi$MGI_ID,mgi_names,result_column="Mouse.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi$human_symbol <- vlookup(mgi$MGI_ID,mgi_names,result_column="Human.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi$human_symbol <- checkGeneSymbols(mgi$human_symbol,unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]

mgi$human_symbol <- vlookup(mgi$MGI_ID,mgi_names,result_column="Human.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi$human_symbol <- checkGeneSymbols(mgi$human_symbol,unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]

mgi <- mgi[-which(is.na(mgi$mouse_symbol)),]  

rm(a,mgi_names,lethal_genes,rows,p)
write.xlsx(mgi,"output/spreadsheets/MGI_genes_with_phenotypes.xlsx")
