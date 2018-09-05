source("Analysis/OMIM_API/API_search_functions.R")

#################searching OMIM for lethal phenotypes######################
search_term<-'(("died as infants") OR ("born died shortly after"~2) OR ("death first year of life"~3) OR ("death infancy"~4) OR ("death infant"~4) OR ("death in utero"~2) OR ("death occurred hours"~4) OR ("death after birth"~3) OR ("delivered dead"~3) OR ("death within birth"~4) OR ("death months life"~4) OR ("died infants"~3) OR ("died days"~3) OR ("died before the age") OR ("died infancy"~3) OR ("died months life"~4) OR  ("dead months life"~4) OR ("death neonatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("death perinatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("lethal before birth"~4) OR ("lethal infantile"~4) OR ("lethal perinatal"~4) OR ("lethal prenatal"~4) OR ("lethal early life"~4) OR ("lethal first year of life"~3) OR ("lethal in utero"~2) OR ("die in utero"~2) OR ("lethal malformation"~3) OR ("lethal neonatal"~3) OR ("no development milestones"~2) OR ("severe lethal"~2) OR ("spontaneous abortion") OR ("the infant died") OR ("fatal infancy"~3) OR ("fetal demise"~3) OR ("embryonic lethal"~3) OR ("early lethal"~2) OR ("die after birth"~2) OR ("Died at weeks"~3) OR ("died at months"~4) OR ("die early childhood"~3) OR ("neonatally lethal"~2) OR ("lethal disorder"~4) OR ("died shortly after birth"~2) OR ("died hours after birth"~5) OR ("neonatal severe"~3) OR (lethal AND congenital) OR (lethal AND prenatal) OR (lethal AND perinatal) OR (lethal AND fetal) OR (lethal AND embryonic) OR (lethal AND neonatal) OR (lethal AND infantile) OR ("embryonically lethal"~3) OR ("prenatally lethal"~3) OR ("infants die"~3) OR ("die within the first year of life") OR ("fetal death") OR ("infants died"~3) OR ("died hours after delivery"~5) OR ("died 1 year"~3) OR ("died at months of age"~3) OR ("died in infancy"~4) OR ("stillborn fetus"~3) OR ("die first year of life"~3) OR ("died at months"~3) OR ("died days"~4) OR ("death at days"~3) OR ("death children"~5) OR ("death occurred month-old"~3) OR ("death early life"~3) OR ("death age months"~4) OR ("death by months"~5) OR ("died hours after birth"~7) OR ("died within months"~2) OR ("infant fatal"~3) OR ("infant died") OR ("early death") OR ("died early age"~3) OR ("died at hours"~2) OR ("died days of life"~3) OR ("death at months"~4))'
search_fields <- c("tx_clinical_features","tx_biochemical_features","tx_description","tx_other_features","tx_genotype_phenotype_correlations")

results<-ldply(lapply(1:length(search_fields),function(x) search_omim(search_term,search_fields[x])),data.frame)


results$omim_field_text<- rep(NA,length(results$mim_numbers))
results$omim_field_text[which(results$search_field=="tx_clinical_features")] <- lapply(results$mim_numbers[which(results$search_field=="tx_clinical_features")],function(x) retrieve_field(x,"clinicalFeatures"))
results$omim_field_text[which(results$search_field=="tx_biochemical_features")] <- lapply(results$mim_numbers[which(results$search_field=="tx_biochemical_features")],function(x) retrieve_field(x,"biochemicalFeatures"))
results$omim_field_text[which(results$search_field=="tx_description")] <- lapply(results$mim_numbers[which(results$search_field=="tx_description")],function(x) retrieve_field(x,"description"))
results$omim_field_text[which(results$search_field=="tx_other_features")] <- lapply(results$mim_numbers[which(results$search_field=="tx_other_features")],function(x) retrieve_field(x,"otherFeatures"))
results$omim_field_text[which(results$search_field=="tx_genotype_phenotype_correlations")] <- lapply(results$mim_numbers[which(results$search_field=="tx_genotype_phenotype_correlations")],function(x) retrieve_field(x,"genotypePhenotypeCorrelations"))

results$mim_entry_type<-vlookup(results$mim_numbers,mapping,result_column = "entry_type",lookup_column="mim_number")

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
results$gene_mim_number<-lapply(results$gene_mim_number,unique)
rm(genemaps,gene_mim_numbers,mapping_keys,inheritances,gene_symbols)



results$gene_mim_number <- lapply(results$gene_mim_number, function(x) unique(strsplit(x,", ")))
results$gene_mim_number <- lapply(results$gene_mim_number,function(x) unique(unlist(x)))

mim_ids<-unique(unlist(results$gene_mim_number[which(results$gene_mim_number!="chromosomal abnormality"&results$gene_mim_number!="genetic basis unknown")]))

length(which(mim_ids%in%universe_df$mim_number))
lethal_genes <- universe_df[which(universe_df$mim_number%in%mim_ids),]

write.xlsx(results,"output/spreadsheets/omim_api_search_hpo_update.xlsx")
save(results, file="output/Data/human_lethal_hits.rda", compress="bzip2")
save(lethal_genes, file="output/Data/human_lethal_genes.rda", compress="bzip2")


#########################writing spreadsheet for Sandra########
load("output/Data/human_lethal_hits.rda")
load("output/Data/human_lethal_genes.rda")

#lethal_genes <- lethal_genes[,-c(1,6,11,16,18,20,22,32,34,36)]
extract_phen_id <- function(phenotype){
  no_phens <- length(phenotype)
  start<- unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-7))
  end<-unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-2))
  return(substr(phenotype,start,end))
}


lethal_genes$omim_phen_id<-lapply(lethal_genes$phenotype,extract_phen_id)

lethal_genes$lethal_phenotype_mim <- lapply(lethal_genes$omim_phen_id,function(x) unique(x[which(x%in%results$mim_numbers)]))

#provisional links, disorders caused by somatic mutations
lethal_genes<-lethal_genes[-which(lengths(lethal_genes$lethal_phenotype)==0),]

#getting info on lethal phenotype from results df, pasting onto lethal_genes df
lethal_genes$lethal_phenotype<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results,lookup_column = "mim_numbers",result_column="preferred_titles"))
lethal_genes$matches<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results,lookup_column = "mim_numbers",result_column="matches"))
lethal_genes$search_field<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results,lookup_column = "mim_numbers",result_column="search_field"))
lethal_genes$omim_field_text<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results,lookup_column = "mim_numbers",result_column="omim_field_text"))
lethal_genes$lethal_phenotype<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results,lookup_column = "mim_numbers",result_column="preferred_titles"))

save(lethal_genes, file="output/Data/human_lethal_genes.rda", compress="bzip2")

#write.xlsx(lethal_genes[,-c(25,26,27)],"output/spreadsheets/human_lethal_genes.xlsx")

###########fixing genes with multiple phenotypes so we know the inheritance of the lETHAL phenotype ###############
load("output/Data/human_lethal_genes.rda")
retrieve_inheritance  <- function(omim_id){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=geneMap&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  inheritance <-unlist(xpathApply(xml, "/omim/entryList/entry/phenotypeMapList/phenotypeMap/phenotypeInheritance", xmlValue))
  return(inheritance)
}
lethal_genes$lethal_inh <- lapply(1:length(lethal_genes$lethal_phenotype_mim),function(x) lapply(lethal_genes$lethal_phenotype_mim[[x]],retrieve_inheritance))
lethal_genes$lethal_inh<-lapply(1:length(lethal_genes$lethal_inh),function(x) lapply(lethal_genes$lethal_inh[[x]],unique))

lethal_genes$AR<-lapply(1:length(lethal_genes$lethal_inh),function(x) 
  ifelse(grepl("Autosomal recessive",lethal_genes$lethal_inh[x]),"Y","N"))
lethal_genes$AD<-lapply(1:length(lethal_genes$lethal_inh),function(x) 
  ifelse(grepl("Autosomal dominant",lethal_genes$lethal_inh[x]),"Y","N"))
lethal_genes$XLd<-lapply(1:length(lethal_genes$lethal_inh),function(x) 
  ifelse(grepl("X-linked dominant",lethal_genes$lethal_inh[x]),"Y","N"))
lethal_genes$XLr<-lapply(1:length(lethal_genes$lethal_inh),function(x) 
  ifelse(grepl("X-linked recessive",lethal_genes$lethal_inh[x]),"Y","N"))
lethal_genes$MT<-lapply(1:length(lethal_genes$lethal_inh),function(x) 
  ifelse(grepl("Mitochondrial",lethal_genes$lethal_inh[x]),"Y","N"))
lethal_genes$XL<-lapply(1:length(lethal_genes$lethal_inh),function(x) 
  ifelse((grepl("X-linked",lethal_genes$lethal_inh[x])&!grepl("recessive",lethal_genes$lethal_inh[x])&!grepl("dominant",lethal_genes$lethal_inh[x])),"Y","N"))


lethal_genes$lethal_inheritance_pattern <- rep(NA,length(lethal_genes$gene))
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$MT=="Y")] <- "MT"
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XLd == "Y")]<- ifelse(is.na(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XLd == "Y")]),"XLd",paste(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XLd == "Y")],"XLd",sep=","))
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XLr == "Y")]<- ifelse(is.na(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XLr == "Y")]),"XLr",paste(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XLr == "Y")],"XLr",sep=","))
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$AR == "Y")]<- ifelse(is.na(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$AR == "Y")]),"AR",paste(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$AR == "Y")],"AR",sep=","))
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$AD == "Y")]<- ifelse(is.na(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$AD == "Y")]),"AD",paste(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$AD == "Y")],"AD",sep=","))
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XL == "Y")]<- ifelse(is.na(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XL == "Y")]),"XL",paste(lethal_genes$lethal_inheritance_pattern[which(lethal_genes$XL == "Y")],"XL",sep=","))

#fix inheritances of phenotype 252010
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$gene=="NDUFA1")]<-"MT,XLd"
lethal_genes$lethal_inheritance_pattern[which(lethal_genes$lethal_inheritance_pattern=="MT,XLd,AR")]<-"MT,AR"


save(lethal_genes, file="output/Data/human_lethal_genes.rda", compress="bzip2")


