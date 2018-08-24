
source("Analysis/OMIM_API/API_search_functions.R")

# more terms 16/7

search_term<-'(("died as infants") OR ("born died shortly after"~2) OR ("death first year of life"~3) OR ("death infancy"~4) OR ("death infant"~4) OR ("death in utero"~2) OR ("death occurred hours"~4) OR ("death after birth"~3) OR ("delivered dead"~3) OR ("death within birth"~4) OR ("death months life"~4) OR ("died infants"~3) OR ("died days"~3) OR ("died before the age") OR ("died infancy"~3) OR ("died months life"~4) OR  ("dead months life"~4) OR ("death neonatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("death perinatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("lethal before birth"~4) OR ("lethal infantile"~4) OR ("lethal perinatal"~4) OR ("lethal prenatal"~4) OR ("lethal early life"~4) OR ("lethal first year of life"~3) OR ("lethal in utero"~2) OR ("die in utero"~2) OR ("lethal malformation"~3) OR ("lethal neonatal"~3) OR ("no development milestones"~2) OR ("severe lethal"~2) OR ("spontaneous abortion") OR ("the infant died") OR ("fatal infancy"~3) OR ("fetal demise"~3) OR ("embryonic lethal"~3) OR ("early lethal"~2) OR ("die after birth"~2) OR ("Died at weeks"~3) OR ("died at months"~4) OR ("die early childhood"~3) OR ("neonatally lethal"~2) OR ("lethal disorder"~4) OR ("died shortly after birth"~2) OR ("died hours after birth"~5) OR ("neonatal severe"~3) OR (lethal AND congenital) OR (lethal AND prenatal) OR (lethal AND perinatal) OR (lethal AND fetal) OR (lethal AND embryonic) OR (lethal AND neonatal) OR (lethal AND infantile) OR ("embryonically lethal"~3) OR ("prenatally lethal"~3) OR ("infants die"~3) OR ("die within the first year of life"))'



my_search   <- paste('search=', search_field,":",search_term, sep='')
my_link  <- 'http://api.omim.org/api/entry/search?'
my_query <- paste(my_link, my_search, "&start=0&limit=10000","&", my_key,sep='')


http://api.omim.org/api/entry/search?search=tx_clinical_features:("died in infancy"~5)&start=0&limit=10000&apiKey=AYosiEyaSfa-CAhAuFIgvg

##############Death before birth############


#Death days after birth
("died hours after birth"~6)

#death weeks/months after birth- Death within the first 24 months of life.- Death in infancy
("died in infancy"~5)
("result in early death")

#death years after birth

#death in early adulthood- Death between the age of 16 and 40 years.

("died at years of age"~6)

("died at age"~3)


("died hours after birth"~6)



#sudden death
("sudden death")










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

