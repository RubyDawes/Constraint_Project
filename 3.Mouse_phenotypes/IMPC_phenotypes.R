
impc<-read.xlsx("Gene_lists/IMPC/IMPC_data.xlsx")
impc$Gene_symbol<- checkGeneSymbols(impc$Gene_symbol,unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
