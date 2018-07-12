library(XML)
my_key <- 'apiKey=AYosiEyaSfa-CAhAuFIgvg'

search_field<-"tx_clinical_features"
search_term<-"+fetal +loss"

my_search   <- paste('search=', search_field,":'",search_term,"'", sep='')
  my_link  <- 'http://api.omim.org/api/entry/search?'
  my_query <- paste(my_link, my_search, "&start=0&limit=10000","&", my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)

mim_numbers<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/mimNumber", xmlValue))
preferred_titles<- unlist(xpathApply(xml, "/omim/searchResponse/entryList/entry/titles/preferredTitle", xmlValue))
matches<- unlist(xpathSApply(xml, "/omim/searchResponse/entryList/entry/matches", xmlValue))

  xml_list<-xmlToList(xml)
  return(xml_list)




hey<-parse_omim(search_omim(search_term,search_field))
