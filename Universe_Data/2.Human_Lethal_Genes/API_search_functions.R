#Functions
my_key <-'apiKey=#####################'
search_omim <- function(search_term,search_field){
  my_search   <- paste('search=', search_field,":",search_term, sep='')
  my_link  <- 'http://api.omim.org/api/entry/search?'
  my_query <- paste(my_link, my_search, "&start=0&limit=10000","&", my_key,sep='')
  print(my_search)
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  mim_numbers<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/mimNumber", xmlValue))
  preferred_titles<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/titles/preferredTitle", xmlValue))
  matches<- unlist(xpathSApply(xml, "/omim/searchResponse/entryList/entry/matches", xmlValue))
  search_results<-data.frame(mim_numbers,preferred_titles,matches,rep(search_field,length(mim_numbers)))
  if (ncol(search_results)==4) {
    names(search_results)<-c("mim_numbers","preferred_titles","matches","search_field")
    return(search_results)
  }
}
retrieve_field <- function(omim_id,search_field){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=text:",search_field,"&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  field <-unlist(xpathApply(xml, "/omim/entryList/entry/textSectionList/textSection/textSectionContent", xmlValue))
  return(field)
}
retrieve_prefix <- function(omim_id){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=geneMap&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  prefix <-unlist(xpathApply(xml, "/omim/entryList/entry/prefix", xmlValue))
  if (length(prefix)==0) {prefix<-NA}
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
search_omim_phrase <- function(search_term,search_field){
  my_search   <- paste('search=', search_term, sep='')
  my_link  <- 'http://api.omim.org/api/entry/search?'
  my_query <- paste(my_link, my_search, "&start=0&limit=10000","&", my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  mim_numbers<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/mimNumber", xmlValue))
  preferred_titles<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/titles/preferredTitle", xmlValue))
  matches<- unlist(xpathSApply(xml, "/omim/searchResponse/entryList/entry/matches", xmlValue))
  search_results<-data.frame(mim_numbers,preferred_titles,matches,rep(search_field,length(mim_numbers)))
  if (ncol(search_results)==4) {
    names(search_results)<-c("mim_numbers","preferred_titles","matches","search_field")
    return(search_results)
  }
}
mapping <- read.csv(file = "Gene_lists/OMIM/mim2gene.txt", sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
names(mapping)<-c("mim_number","entry_type", "entrez_gene_id", "hgnc_symbol")

retrieve_inheritance  <- function(omim_id){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=geneMap&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  inheritance <-unlist(xpathApply(xml, "/omim/entryList/entry/phenotypeMapList/phenotypeMap/phenotypeInheritance", xmlValue))
  return(inheritance)
}