length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z<3.09)) #202 AD genes have LoF constraint but no missense constraint- pool 1
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z>=3.09)) #60 AD genes have missense constraint but no LoF constraint- pool 2
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z>=3.09)) #182 AD genes have missense constraint AND LoF constraint- pool 3
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z<3.09)) #387 AD genes have neither missense constraint nor LoF constraint- pool 4

my_key <- 'apiKey=AYosiEyaSfa-CAhAuFIgvg'


retrieve_variants <- function(omim_id){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=allelicVariantList&include=existFlags&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  phenotypes<- unlist(xpathApply(xml, "/omim/entryList/entry/allelicVariantList/allelicVariant/name", xmlValue))
  if (length(phenotypes)==0){
    phenotypes <- "no allelic variants"
    mutations <- "no allelic variants"
    text <- "no allelic variants"
    search_results<- data.frame(phenotypes,mutations,text)
    names(search_results)<- c("phenotypes","mutations","text")
    return(search_results)}
  if (length(which(phenotypes=="REMOVED FROM DATABASE"))>0) {
    phenotypes <- phenotypes[-which(phenotypes=="REMOVED FROM DATABASE")]
  }
  if (length(which(grepl("MOVED TO",phenotypes)))>0){
    phenotypes <- phenotypes[-which(grepl("MOVED TO",phenotypes))]
  }
  mutations<-unlist(xpathApply(xml, "/omim/entryList/entry/allelicVariantList/allelicVariant/mutations", xmlValue))
  missense_mutations <- unlist(lapply(1:length(mutations),function(x) grepl("Substitution",text[x],ignore.case = TRUE)))
  lof_mutations <- unlist(lapply(1:length(mutations),function(x) grepl("Frameshift",text[x],ignore.case = TRUE)|grepl("premature termination codon",text[x],ignore.case=TRUE)|grepl("nonsense",text[x],ignore.case=TRUE)|grepl("-to-ter",text[x],ignore.case=TRUE)))
  text <-unlist(xpathApply(xml, "/omim/entryList/entry/allelicVariantList/allelicVariant/text", xmlValue))
  if (length(which(grepl("MOLECULAR DEFECT UNKNOWN",phenotypes)))>0){
    a<-which(grepl("MOLECULAR DEFECT UNKNOWN",phenotypes))
    phenotypes <- phenotypes[-a]
    text <- text[-a]
  }
  search_results<- data.frame(phenotypes,mutations,missense_mutations,lof_mutations,text)
  names(search_results)<- c("phenotypes","mutations","missense_mutations","lof_mutations","text")
  return(search_results)
}
###POOL1 AD genes have LoF constraint but no missense constraint- pool 1
#fix later
universe_df$mim_number[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z<3.09)][99]<-"600049"

pool1<-data.frame(mim_ids <-universe_df$mim_number[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z<3.09)], genes <- universe_df$gene[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z<3.09)])
names(pool1)<-c("mim_ids","genes")

pool1_results <- lapply(pool1$mim_ids,retrieve_variants)

pool1$allelic_variants <- unlist(lapply(1:length(pool1$mim_ids),function(x) length(pool1_results[[x]]$text)))
pool1$substitutions <- unlist(lapply(1:length(pool1$mim_ids),function(x) length(which(grepl("Substitution",pool1_results[[x]]$text,ignore.case=TRUE)|grepl("missense",pool1_results[[x]]$text,ignore.case=TRUE)))))
pool1$frameshifts <- unlist(lapply(1:length(pool1$mim_ids),function(x) length(which(grepl("Frameshift",pool1_results[[x]]$text,ignore.case=TRUE)|grepl("premature termination codon",pool1_results[[x]]$text,ignore.case=TRUE)|grepl("nonsense",pool1_results[[x]]$text,ignore.case=TRUE)))))
pool1$unsure <- lapply(1:length(pool1$mim_ids),function(x) pool1$allelic_variants[x]-pool1$substitutions[x]-pool1$frameshifts[x])

sum(pool1$substitutions)/sum(pool1$allelic_variants)
sum(pool1$frameshifts)/sum(pool1$allelic_variants)

###POOL2 AD genes have missense constraint but no LoF constraint- pool 2
pool2<-data.frame(mim_ids <-universe_df$mim_number[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z>=3.09)], genes <- universe_df$gene[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z>=3.09)])
names(pool2)<-c("mim_ids","genes")

pool2_results <- lapply(pool2$mim_ids,retrieve_variants)

pool2$allelic_variants <- unlist(lapply(1:length(pool2$mim_ids),function(x) length(pool2_results[[x]]$text)))
pool2$substitutions <- unlist(lapply(1:length(pool2$mim_ids),function(x) length(which(grepl("Substitution",pool2_results[[x]]$text,ignore.case=TRUE)|grepl("missense",pool2_results[[x]]$text,ignore.case=TRUE)))))
pool2$frameshifts <- unlist(lapply(1:length(pool2$mim_ids),function(x) length(which(grepl("Frameshift",pool2_results[[x]]$text,ignore.case=TRUE)|grepl("premature termination codon",pool2_results[[x]]$text,ignore.case=TRUE)|grepl("nonsense",pool2_results[[x]]$text,ignore.case=TRUE)))))
pool2$unsure <- lapply(1:length(pool2$mim_ids),function(x) pool2$allelic_variants[x]-pool2$substitutions[x]-pool2$frameshifts[x])

sum(pool2$frameshifts)/sum(pool2$allelic_variants)
sum(pool2$substitutions)/sum(pool2$allelic_variants)

###POOL3 182 AD genes have missense constraint AND LoF constraint- pool 3
pool3<-data.frame(mim_ids <-universe_df$mim_number[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z>=3.09)], genes <- universe_df$gene[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z>=3.09)])
names(pool3)<-c("mim_ids","genes")

pool3_results <- lapply(pool3$mim_ids,retrieve_variants)

pool3$allelic_variants <- unlist(lapply(1:length(pool3$mim_ids),function(x) length(pool3_results[[x]]$text)))
pool3$substitutions <- unlist(lapply(1:length(pool3$mim_ids),function(x) length(which(grepl("Substitution",pool3_results[[x]]$text,ignore.case=TRUE)|grepl("missense",pool3_results[[x]]$text,ignore.case=TRUE)))))
pool3$frameshifts <- unlist(lapply(1:length(pool3$mim_ids),function(x) length(which(grepl("Frameshift",pool3_results[[x]]$text,ignore.case=TRUE)|grepl("premature termination codon",pool3_results[[x]]$text,ignore.case=TRUE)|grepl("nonsense",pool3_results[[x]]$text,ignore.case=TRUE)))))
pool3$unsure <- lapply(1:length(pool3$mim_ids),function(x) pool3$allelic_variants[x]-pool3$substitutions[x]-pool3$frameshifts[x])

sum(pool3$frameshifts)/sum(pool3$allelic_variants)
sum(pool3$substitutions)/sum(pool3$allelic_variants)

###POOL4 387 AD genes have neither missense constraint nor LoF constraint- pool 4
pool4<-data.frame(mim_ids <-universe_df$mim_number[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z<3.09)], genes <- universe_df$gene[which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z<3.09)])
names(pool4)<-c("mim_ids","genes")

#fix later
pool4<-pool4[-144,]


pool4_results <- lapply(pool4$mim_ids,retrieve_variants)
for (i in 1:length(pool4$mim_ids)){
  retrieve_variants(pool4$mim_ids[i])
}
pool4$allelic_variants <- unlist(lapply(1:length(pool4$mim_ids),function(x) length(pool4_results[[x]]$text)))
pool4$substitutions <- unlist(lapply(1:length(pool4$mim_ids),function(x) length(which(grepl("Substitution",pool4_results[[x]]$text,ignore.case=TRUE)|grepl("missense",pool4_results[[x]]$text,ignore.case=TRUE)))))
pool4$frameshifts <- unlist(lapply(1:length(pool4$mim_ids),function(x) length(which(grepl("Frameshift",pool4_results[[x]]$text,ignore.case=TRUE)|grepl("premature termination codon",pool4_results[[x]]$text,ignore.case=TRUE)|grepl("nonsense",pool4_results[[x]]$text,ignore.case=TRUE)))))
pool4$unsure <- lapply(1:length(pool4$mim_ids),function(x) pool4$allelic_variants[x]-pool4$substitutions[x]-pool4$frameshifts[x])

sum(pool4$frameshifts)/sum(pool4$allelic_variants)
sum(pool4$substitutions)/sum(pool4$allelic_variants)


#202 AD genes have LoF constraint but no missense constraint- pool 1
length(which(pool1$substitutions>0))/length(pool1$mim_ids)
length(which(pool1$frameshifts>0))/length(pool1$mim_ids)
#60 AD genes have missense constraint but no LoF constraint- pool 2
length(which(pool2$substitutions>0))/length(pool2$mim_ids)
length(which(pool2$frameshifts>0))/length(pool2$mim_ids)
#182 AD genes have missense constraint AND LoF constraint- pool 3
length(which(pool3$substitutions>0))/length(pool3$mim_ids)
length(which(pool3$frameshifts>0))/length(pool3$mim_ids)
#387 AD genes have neither missense constraint nor LoF constraint- pool 4
length(which(pool4$substitutions>0))/length(pool4$mim_ids)
length(which(pool4$frameshifts>0))/length(pool4$mim_ids)

library("hpoPlot")
hp <- get.ontology("Gene_lists/HPO/hp.obo",qualifier="HP")
universe_df$hpo_terms[which(universe_df$mim_number=="605452")]
hey<-c("HP:0000992","HP:0004322","HP:0004802","HP:0001923","HP:0000006","HP:0001480","HP:0002153","HP:0001034")
get.ancestors(hp,hey)

general_cats <- c("HP:0000707","HP:0000478","HP:0000152","HP:0000119","HP:0000924","HP:0001939","HP:0003011",
                  "HP:0001871","HP:0001626","HP:0002664","HP:0002715","HP:0001574","HP:0040064","HP:0025031",
                  "HP:0000598","HP:0000818","HP:0002086","HP:0001197","HP:0003549","HP:0000769","HP:0001507","HP:0045027")
general_cats_names <- get.shortened.names(hp,general_cats)

cell_ess_hpo_anc <- lapply(cell_ess_hpo,function(x) get.ancestors(hp,x))
library(httr)

clinvar_txt <- read.csv("Gene_lists/ClinVar/variant_summary.txt",sep="\t",stringsAsFactors = FALSE)
  
retrieve_variants_clinvar<- function(omim_id){
  my_mim   <- paste('mimNumber=', omim_id, sep='')
  my_link  <- 'http://api.omim.org/api/entry?'
  my_query <- paste(my_link, my_mim, "&include=allelicVariantList&include=existFlags&",my_key,sep='')
  xml<-xmlTreeParse(my_query, useInternalNodes=TRUE)
  clinvar <- unlist(xpathApply(xml, "/omim/entryList/entry/allelicVariantList/allelicVariant/clinvarAccessions",xmlValue))
  clinvar<-unlist(strsplit(clinvar,";;;"))
  clinvar_query<-unlist(lapply(1:length(clinvar),function(x) paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=clinvarset&id=",clinvar[x],sep="")))
  clinvar_xml <- unlist(lapply(1:length(clinvar_query),function(x) xmlTreeParse(rawToChar(GET(clinvar_query[x])$content),useInternalNodes=TRUE)))
  clinvar_consequence<-lapply(1:length(clinvar_xml),function(x) xpathApply(clinvar_xml[[x]],"/ClinVarResult-Set/ClinVarSet/ReferenceClinVarAssertion/MeasureSet[@Type='Variant']/Measure/AttributeSet/Attribute[@Type='MolecularConsequence']"))
  clinvar_consequence<-lapply(1:length(clinvar_consequence),function(x) unlist(lapply(clinvar_consequence[[x]],xmlValue),recursive = FALSE))

  mim_phen_id <- unlist(lapply(1:length(clinvar),function(x) vlookup(clinvar[x],clinvar_txt,result_column="PhenotypeIDS",lookup_column="RCVaccession")))
  mim_phen <- unlist(lapply(1:length(clinvar),function(x) vlookup(clinvar[x],clinvar_txt,result_column="PhenotypeList",lookup_column="RCVaccession")))
  results<- data.frame(clinvar,clinvar_consequence,mim_phen_id,mim_phen)
  names(results)<- c("clinvar","clinvar_consequence","mim_phen_id","mim_phen")
  return(results)
}
xpathApply(clinvar_xml[[1]],'/ClinVarResult-Set/ClinVarSet/ReferenceClinVarAssertion/TraitSet[@Type="Disease"]/Trait/Name/ElementValue[@Type="Preferred"]',xmlValue)
#clinvar stuff
pool1_clinvar <- unlist(lapply(pool1$mim_ids,retrieve_variants_clinvar))
for (i in 1:length(pool1$mim_ids)){
  retrieve_variants_clinvar(pool1$mim_ids[i])
}


  