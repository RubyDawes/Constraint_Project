#CREDIT: functions adapted from "https://davetang.org/muse/2015/03/17/getting-started-with-the-omim-api/"
#load library

source("OMIM_API/API_search_functions.R")

################searching omim for entries related to foetal/perinatal lethality#################
search_term<-"((lethal AND congenital) OR (lethal AND prenatal) OR (lethal AND perinatal) OR (lethal AND fetal) OR (lethal AND embryonic) OR (lethal AND neonatal) OR (lethal AND infantile))"
search_fields <- c("tx_clinical_features","tx_biochemical_features","tx_description",
                   "tx_other_features","tx_genotype_phenotype_correlations")
results<- ldply(lapply(search_fields, function(y) search_omim(search_term,y)),data.frame)

results$mim_entry_type<-vlookup(results$mim_numbers,mapping,result_column = "entry_type",lookup_column="mim_number")

results$omim_field_text<- rep(NA,length(results$mim_numbers))
results$omim_field_text[which(results$search_field=="tx_clinical_features")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_clinical_features")],function(x) retrieve_field(x,"clinicalFeatures")))
results$omim_field_text[which(results$search_field=="tx_biochemical_features")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_biochemical_features")],function(x) retrieve_field(x,"biochemicalFeatures")))
results$omim_field_text[which(results$search_field=="tx_description")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_description")],function(x) retrieve_field(x,"description")))
results$omim_field_text[which(results$search_field=="tx_other_features")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_other_features")],function(x) retrieve_field(x,"otherFeatures")))
results$omim_field_text[which(results$search_field=="tx_genotype_phenotype_correlations")] <- unlist(lapply(results$mim_numbers[which(results$search_field=="tx_genotype_phenotype_correlations")],function(x) retrieve_field(x,"genotypePhenotypeCorrelations")))

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

results1<- read.xlsx("output/spreadsheets/curated_lethal_list.xlsx")
mim_ids<-unique(unlist(strsplit(gsub(" ","",paste(results$gene_mim_number[which(results$gene_mim_number!="chromosomal abnormality"&results$gene_mim_number!="genetic basis unknown")],collapse=","),fixed=TRUE),",")))
length(which(mim_ids%in%universe_df$mim_number))

lethal_genes <- universe_df[which(universe_df$mim_number%in%mim_ids),]

#expanding existing search to new fields
search_term<-"((lethal AND congenital) OR (lethal AND prenatal) OR (lethal AND perinatal) OR (lethal AND fetal) OR (lethal AND embryonic) OR (lethal AND neonatal) OR (lethal AND infantile) OR (stillbirth) )"
search_fields<- c("tx_clinical_management","tx_cytogenetics","tx_diagnosis","tx_genetic_variability","tx_genotype","tx_heterogeneity","tx_inheritance","tx_molecular_genetics","tx_pathogenesis","tx_phenotype")
results2<- ldply(lapply(search_fields, function(y) search_omim(search_term,y)),data.frame)
results2_unique<-results2[which(!results2$mim_numbers%in%results1$mim_numbers),]
results2_unique$omim_field_text<- rep(NA,length(results2_unique$mim_numbers))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_clinical_management")] <- 
  lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_clinical_management")],function(x) retrieve_field(x,"clinicalManagement"))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_cytogenetics")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_cytogenetics")],function(x) retrieve_field(x,"cytogenetics")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_diagnosis")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_diagnosis")],function(x) retrieve_field(x,"diagnosis")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_genetic_variability")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_genetic_variability")],function(x) retrieve_field(x,"geneticVariability")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_genotype")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_genotype")],function(x) retrieve_field(x,"genotype")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_heterogeneity")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_heterogeneity")],function(x) retrieve_field(x,"heterogeneity")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_inheritance")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_inheritance")],function(x) retrieve_field(x,"inheritance")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_molecular_genetics")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_molecular_genetics")],function(x) retrieve_field(x,"molecularGenetics")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_pathogenesis")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_pathogenesis")],function(x) retrieve_field(x,"pathogenesis")))
results2_unique$omim_field_text[which(results2_unique$search_field=="tx_phenotype")] <- 
  unlist(lapply(results2_unique$mim_numbers[which(results2_unique$search_field=="tx_phenotype")],function(x) retrieve_field(x,"phenotype")))

results2_unique$prefix<-unlist(lapply(results2_unique$mim_numbers, retrieve_prefix))
genemaps<-  lapply(results2_unique$mim_numbers[which(results2_unique$prefix=="#")],retrieve_genemap)
gene_mim_numbers <- c()
mapping_keys <- c()
inheritances <- c()
gene_symbols <- c()
for (i in 1:length(which(results2_unique$prefix=="#"))) {
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

results2_unique$gene_mim_number[which(results2_unique$prefix=="#")]<- gene_mim_numbers
results2_unique$phenotype_mapping_keys[which(results2_unique$prefix=="#")]<- mapping_keys
results2_unique$inheritance[which(results2_unique$prefix=="#")]<- inheritances
results2_unique$gene_symbol[which(results2_unique$prefix=="#")]<- gene_symbols
results2_unique$gene_mim_number[which(results2_unique$prefix=="*")]<-as.character(results2_unique$mim_numbers[which(results2_unique$prefix=="*")])
results2_unique$gene_mim_number[which(results2_unique$prefix=="%"|results2_unique$mim_entry_type=="predominantly phenotypes")]<-"genetic basis unknown"
results2_unique$gene_mim_number[which(is.na(results2_unique$gene_mim_number))]<-"chromosomal abnormality"
rm(genemaps,gene_mim_numbers,mapping_keys,inheritances,gene_symbols)

mim_ids<-unique(unlist(strsplit(gsub(" ","",paste(results2_unique$gene_mim_number[which(results2_unique$gene_mim_number!="chromosomal abnormality"&results2_unique$gene_mim_number!="genetic basis unknown")],collapse=","),fixed=TRUE),",")))
length(which(mim_ids%in%universe_df$mim_number))
lethal_genes <- universe_df[which(universe_df$mim_number%in%mim_ids),]


# more terms 16/7

search_term<-'(("died as infants") OR ("born died shortly after"~2) OR ("death first year of life"~3) OR ("death infancy"~4) OR ("death infant"~4) OR ("death in utero"~2) OR ("death occurred hours"~4) OR ("death after birth"~3) OR ("delivered dead"~3) OR ("death within birth"~4) OR ("death months life"~4) OR ("died infants"~3) OR ("died days"~3) OR ("died before the age") OR ("died infancy"~3) OR ("died months life"~4) OR  ("dead months life"~4) OR ("death neonatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("death perinatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("lethal before birth"~4) OR ("lethal infantile"~4) OR ("lethal perinatal"~4) OR ("lethal prenatal"~4) OR ("lethal early life"~4) OR ("lethal first year of life"~3) OR ("lethal in utero"~2) OR ("die in utero"~2) OR ("lethal malformation"~3) OR ("lethal neonatal"~3) OR ("no development milestones"~2) OR ("severe lethal"~2) OR ("spontaneous abortion") OR ("the infant died") OR ("fatal infancy"~3) OR ("fetal demise"~3) OR ("embryonic lethal"~3) OR ("early lethal"~2) OR ("die after birth"~2) OR ("Died at weeks"~3) OR ("died at months"~4) OR ("die early childhood"~3) OR ("neonatally lethal"~2) OR ("lethal disorder"~4) OR ("died shortly after birth"~2) OR ("died hours after birth"~5) OR ("neonatal severe"~3) OR (lethal AND congenital) OR (lethal AND prenatal) OR (lethal AND perinatal) OR (lethal AND fetal) OR (lethal AND embryonic) OR (lethal AND neonatal) OR (lethal AND infantile) OR ("embryonically lethal"~3) OR ("prenatally lethal"~3) OR ("infants die"~3) OR ("die within the first year of life"))'
search_fields <- c("tx_clinical_features","tx_biochemical_features","tx_description","tx_other_features","tx_genotype_phenotype_correlations")
results3<-ldply(lapply(1:length(search_fields),function(x) search_omim(search_term,search_fields[x])),data.frame)


results3$omim_field_text<- rep(NA,length(results3$mim_numbers))
results3$omim_field_text[which(results3$search_field=="tx_clinical_features")] <- lapply(results3$mim_numbers[which(results3$search_field=="tx_clinical_features")],function(x) retrieve_field(x,"clinicalFeatures"))
results3$omim_field_text[which(results3$search_field=="tx_biochemical_features")] <- lapply(results3$mim_numbers[which(results3$search_field=="tx_biochemical_features")],function(x) retrieve_field(x,"biochemicalFeatures"))
results3$omim_field_text[which(results3$search_field=="tx_description")] <- lapply(results3$mim_numbers[which(results3$search_field=="tx_description")],function(x) retrieve_field(x,"description"))
results3$omim_field_text[which(results3$search_field=="tx_other_features")] <- lapply(results3$mim_numbers[which(results3$search_field=="tx_other_features")],function(x) retrieve_field(x,"otherFeatures"))
results3$omim_field_text[which(results3$search_field=="tx_genotype_phenotype_correlations")] <- lapply(results3$mim_numbers[which(results3$search_field=="tx_genotype_phenotype_correlations")],function(x) retrieve_field(x,"genotypePhenotypeCorrelations"))

results3$mim_entry_type<-vlookup(results3$mim_numbers,mapping,result_column = "entry_type",lookup_column="mim_number")

results3$prefix<-unlist(lapply(results3$mim_numbers, retrieve_prefix))
genemaps<-  lapply(results3$mim_numbers[which(results3$prefix=="#")],retrieve_genemap)

gene_mim_numbers <- c()
mapping_keys <- c()
inheritances <- c()
gene_symbols <- c()
for (i in 1:length(which(results3$prefix=="#"))) {
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

results3$gene_mim_number[which(results3$prefix=="#")]<- gene_mim_numbers
results3$phenotype_mapping_keys[which(results3$prefix=="#")]<- mapping_keys
results3$inheritance[which(results3$prefix=="#")]<- inheritances
results3$gene_symbol[which(results3$prefix=="#")]<- gene_symbols
results3$gene_mim_number[which(results3$prefix=="*")]<-as.character(results3$mim_numbers[which(results3$prefix=="*")])
results3$gene_mim_number[which(results3$prefix=="%"|results3$mim_entry_type=="predominantly phenotypes")]<-"genetic basis unknown"
results3$gene_mim_number[which(is.na(results3$gene_mim_number))]<-"chromosomal abnormality"
results3$gene_mim_number<-lapply(results3$gene_mim_number,unique)
rm(genemaps,gene_mim_numbers,mapping_keys,inheritances,gene_symbols)



results3$gene_mim_number <- lapply(results3$gene_mim_number, function(x) unique(strsplit(x,", ")))
results3$gene_mim_number <- lapply(results3$gene_mim_number,function(x) unique(unlist(x)))

results3$mouse_phen <-lapply(results3$gene_mim_number,function(x) na.omit(vlookup(x,lethal_genes,result_column = "lethal_mouse",lookup_column="mim_number")))

exclude<- read.xlsx("output/spreadsheets/omim_api_exclude.xlsx")
results3<-results3[-which(results3$mim_numbers%in%exclude$mim_numbers[1:20]),]

mim_ids<-unique(unlist(results3$gene_mim_number[which(results3$gene_mim_number!="chromosomal abnormality"&results3$gene_mim_number!="genetic basis unknown")]))

length(which(mim_ids%in%universe_df$mim_number))
lethal_genes <- universe_df[which(universe_df$mim_number%in%mim_ids),]

write.xlsx(results3,"output/spreadsheets/omim_api_search_results_better.xlsx")
save(results3, file="output/Data/human_lethal_hits.rda", compress="bzip2")
save(lethal_genes, file="output/Data/human_lethal_genes.rda", compress="bzip2")


#########################writing spreadsheet for Sandra########
load("output/Data/human_lethal_hits.rda")
load("output/Data/human_lethal_genes.rda")

lethal_genes <- lethal_genes[,-c(1,6,11,16,18,20,22,32,34,36)]
extract_phen_id <- function(phenotype){
  no_phens <- length(phenotype)
  start<- unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-7))
  end<-unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-2))
  return(substr(phenotype,start,end))
}

extract_phen_id <- function(phenotype){
  no_phens <- length(phenotype)
  start<- unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-7))
  end<-unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-2))
  return(substr(phenotype,start,end))
}

lethal_genes$omim_phen_id<-lapply(lethal_genes$phenotype,extract_phen_id)

lethal_genes$lethal_phenotype_mim <- lapply(lethal_genes$omim_phen_id,function(x) unique(x[which(x%in%results3$mim_numbers)]))

#provisional links, disorders caused by somatic mutations
lethal_genes<-lethal_genes[-which(lengths(lethal_genes$lethal_phenotype)==0),]

#getting info on lethal phenotype from results df, pasting onto lethal_genes df
lethal_genes$lethal_phenotype<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results3,lookup_column = "mim_numbers",result_column="preferred_titles"))
lethal_genes$matches<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results3,lookup_column = "mim_numbers",result_column="matches"))
lethal_genes$search_field<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results3,lookup_column = "mim_numbers",result_column="search_field"))
lethal_genes$omim_field_text<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results3,lookup_column = "mim_numbers",result_column="omim_field_text"))
lethal_genes$lethal_phenotype<-lapply(lethal_genes$lethal_phenotype_mim,function(x) 
  vlookup(x,results3,lookup_column = "mim_numbers",result_column="preferred_titles"))

save(lethal_genes, file="output/Data/human_lethal_genes.rda", compress="bzip2")

#write.xlsx(lethal_genes[,-c(40:45)],"output/spreadsheets/human_lethal_genes.xlsx")

