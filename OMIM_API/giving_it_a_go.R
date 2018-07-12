#CREDIT: functions adapted from "https://davetang.org/muse/2015/03/17/getting-started-with-the-omim-api/"
#load library


search_omim <- function(search_term,search_field){
  my_search   <- paste('search=', search_field,":",search_term, sep='')
  my_link  <- 'http://api.omim.org/api/entry/search?'
  my_query <- paste(my_link, my_search, "&start=0&limit=10000","&", my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  mim_numbers<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/mimNumber", xmlValue))
  preferred_titles<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/titles/preferredTitle", xmlValue))
  matches<- unlist(xpathSApply(xml, "/omim/searchResponse/entryList/entry/matches", xmlValue))
  search_results<-data.frame(mim_numbers,preferred_titles,matches,rep(search_field,length(mim_numbers)))
  names(search_results)<-c("mim_numbers","preferred_titles","matches","search_field")
  return(search_results)
}

#searching omim for entries related to foetal/perinatal lethality
search_term<-"((lethal AND congenital) OR (lethal AND prenatal) OR (lethal AND perinatal) OR (lethal AND fetal) OR (lethal AND embryonic) OR (lethal AND neonatal) OR (lethal AND infantile) OR (stillbirth) )"
search_fields <- c("tx_clinical_features","tx_biochemical_features","tx_description",
                   "tx_other_features","tx_genotype_phenotype_correlations")
results<- ldply(lapply(search_fields, function(y) search_omim(search_term,y)),data.frame)

mapping <- read.csv(file = "Gene_lists/OMIM/mim2gene.txt", sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
names(mapping)<-c("mim_number","entry_type", "entrez_gene_id", "hgnc_symbol")
results$mim_entry_type<-vlookup(results$mim_numbers,mapping,result_column = "entry_type",lookup_column="mim_number")



include=text:clinicalFeatures

retrieve_field <- function(omim_id,search_field){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=text:",search_field,"&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  field <-unlist(xpathApply(xml, "/omim/entryList/entry/textSectionList/textSection/textSectionContent", xmlValue))
  return(field)
}

results$omim_field_text<- rep(NA,length(results$mim_numbers))
results$omim_field_text[which(results$search_field=="tx_clinical_features")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_clinical_features")],function(x) retrieve_field(x,"clinicalFeatures")))
results$omim_field_text[which(results$search_field=="tx_biochemical_features")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_biochemical_features")],function(x) retrieve_field(x,"biochemicalFeatures")))
results$omim_field_text[which(results$search_field=="tx_description")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_description")],function(x) retrieve_field(x,"description")))
results$omim_field_text[which(results$search_field=="tx_other_features")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_other_features")],function(x) retrieve_field(x,"otherFeatures")))
results$omim_field_text[which(results$search_field=="tx_genotype_phenotype_correlations")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_genotype_phenotype_correlations")],function(x) retrieve_field(x,"genotypePhenotypeCorrelations")))


retrieve_prefix <- function(omim_id){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=geneMap&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  prefix <-unlist(xpathApply(xml, "/omim/entryList/entry/prefix", xmlValue))
  return(prefix)}

retrieve_genemap <- function(omim_id){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=geneMap&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  gene_mim_numbers<- unlist(xpathApply(xml, "/omim/entryList/entry/phenotypeMapList/phenotypeMap/mimNumber", xmlValue))
  phenotype_mapping_key <-unlist(xpathApply(xml, "/omim/entryList/entry/phenotypeMapList/phenotypeMap/phenotypeMappingKey", xmlValue))
  inheritance <-unlist(xpathApply(xml, "/omim/entryList/entry/phenotypeMapList/phenotypeMap/phenotypeInheritance", xmlValue))
  gene_symbol <-unlist(xpathApply(xml, "/omim/entryList/entry/phenotypeMapList/phenotypeMap/geneSymbols", xmlValue))
  search_results<- data.frame(gene_mim_numbers,phenotype_mapping_key,inheritance,gene_symbol)
  if (ncol(search_results)==4){
    names(search_results)<-c("gene_mim_numbers","phenotype_mapping_key","inheritance","gene_symbol")
    return(search_results)
  }
  else {return("no associated gene")}
  }


results$prefix<-unlist(lapply(results$mim_numbers, retrieve_prefix))
genemaps<-  lapply(results$mim_numbers[which(results$prefix=="#")],retrieve_genemap)
gene_mim_numbers <- c()
mapping_keys <- c()
inheritances <- c()
gene_symbols <- c()
for (i in 1:length(which(results$prefix=="#"))) {
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

results$gene_mim_number[which(results$prefix=="#")]<- gene_mim_numbers
results$phenotype_mapping_keys[which(results$prefix=="#")]<- mapping_keys
results$inheritance[which(results$prefix=="#")]<- inheritances
results$gene_symbol[which(results$prefix=="#")]<- gene_symbols
results$gene_mim_number[which(results$prefix=="*")]<-as.character(results$mim_numbers[which(results$prefix=="*")])
results$gene_mim_number[which(results$prefix=="%"|results$mim_entry_type=="predominantly phenotypes")]<-"genetic basis unknown"
results$gene_mim_number[which(is.na(results$gene_mim_number))]<-"chromosomal abnormality"
rm(genemaps,gene_mim_numbers,mapping_keys,inheritances,gene_symbols)

write.xlsx(results,"output/spreadsheets/omim_api_search.xlsx")

mim_ids<-unique(unlist(strsplit(gsub(" ","",paste(results$gene_mim_number[which(results$gene_mim_number!="chromosomal abnormality"&results$gene_mim_number!="genetic basis unknown")],collapse=","),fixed=TRUE),",")))
length(which(mim_ids%in%universe_df$mim_number))



lethal_genes <- universe_df[which(universe_df$mim_number%in%mim_ids),]


