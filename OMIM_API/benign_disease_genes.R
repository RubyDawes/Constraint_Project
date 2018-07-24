source("OMIM_API/API_search_functions.R")

benign_search_fields<-c("title","tx_description")
benign_hits<-ldply(lapply(benign_search_fields,function(x) search_omim("benign",x)),data.frame)


benign_hits$mim_entry_type<-vlookup(benign_hits$mim_numbers,mapping,result_column = "entry_type",lookup_column="mim_number")

benign_hits$omim_field_text<- rep(NA,length(benign_hits$mim_numbers))
benign_hits$omim_field_text[which(benign_hits$search_field=="tx_description")] <- unlist(lapply(benign_hits$mim_numbers[which(benign_hits$search_field=="tx_description")],function(x) retrieve_field(x,"description")))

benign_hits$prefix<-unlist(lapply(benign_hits$mim_numbers, retrieve_prefix))
genemaps<-  lapply(benign_hits$mim_numbers[which(benign_hits$prefix=="#")],retrieve_genemap)
gene_mim_numbers <- c()
mapping_keys <- c()
inheritances <- c()
gene_symbols <- c()
for (i in 1:length(which(benign_hits$prefix=="#"))) {
  if (typeof(genemaps[[i]])=="character"){
    gene_mim_numbers[i] <- NA
    mapping_keys[i]<-NA
    inheritances[i]<-NA
    gene_symbols[i]<-genemaps[[i]]
  }
  else {
    gene_mim_numbers[i]<-paste(as.character(genemaps[[i]][1]$gene_mim_numbers),collapse=", ")
    mapping_keys[i]<-paste(as.character(genemaps[[i]][2]$phenotype_mapping_key),collapse=", ")
    inheritances[i]<-paste(as.character(genemaps[[i]][3]$inheritance),collapse=", ")
    gene_symbols[i]<-paste(as.character(genemaps[[i]][4]$gene_symbol),collapse=", ")
  }}
rm(i)

benign_hits$gene_mim_number[which(benign_hits$prefix=="#")]<- gene_mim_numbers
benign_hits$phenotype_mapping_keys[which(benign_hits$prefix=="#")]<- mapping_keys
benign_hits$inheritance[which(benign_hits$prefix=="#")]<- inheritances
benign_hits$gene_symbol[which(benign_hits$prefix=="#")]<- gene_symbols
benign_hits$gene_mim_number[which(benign_hits$prefix=="*")]<-as.character(benign_hits$mim_numbers[which(benign_hits$prefix=="*")])
benign_hits$gene_mim_number[which(benign_hits$prefix=="%"|benign_hits$mim_entry_type=="predominantly phenotypes")]<-"genetic basis unknown"
benign_hits$gene_mim_number[which(is.na(benign_hits$gene_mim_number))]<-"chromosomal abnormality"
rm(genemaps,gene_mim_numbers,mapping_keys,inheritances,gene_symbols)

mim_ids<-unique(unlist(strsplit(gsub(" ","",paste(benign_hits$gene_mim_number[which(benign_hits$gene_mim_number!="chromosomal abnormality"&benign_hits$gene_mim_number!="genetic basis unknown")],collapse=","),fixed=TRUE),",")))
length(which(mim_ids%in%universe_df$mim_number))

benign_genes <- universe_df[which(universe_df$mim_number%in%mim_ids),]
benign_only <- benign_genes[which(lengths(benign_genes$phenotype)==1),]

