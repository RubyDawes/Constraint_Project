

extract_phen_id <- function(phenotype){
  no_phens <- length(phenotype)
  start<- unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-7))
  end<-unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-2))
  return(substr(phenotype,start,end))
}

universe_df$omim_phen_id<-lapply(universe_df$phenotype,extract_phen_id)

clinvar_txt <- read.csv("Gene_lists/ClinVar/variant_summary.txt",sep="\t",stringsAsFactors = FALSE)
clinvar_txt<- clinvar_txt[which(clinvar_txt$ClinSigSimple==1&clinvar_txt$Assembly=="GRCh38"),]

retrieve_variant_consequence<- function(RCVaccession){
  clinvar_query<-paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=clinvarset&id=",RCVaccession,sep="")
  clinvar_xml <- xmlTreeParse(rawToChar(GET(clinvar_query)$content),useInternalNodes=TRUE)
  clinvar_consequence<-xpathApply(clinvar_xml,"/ClinVarResult-Set/ClinVarSet/ReferenceClinVarAssertion/MeasureSet[@Type='Variant']/Measure/AttributeSet/Attribute[@Type='MolecularConsequence']")
  clinvar_consequence<-unlist(lapply(clinvar_consequence,xmlValue),recursive = FALSE)
  return(clinvar_consequence)
}
find_phenotype_causal_variants <- function(phen_mim_no){
  mim_no <- paste("OMIM:",phen_mim_no,sep="")
  clinvar_accessions <- clinvar_txt$RCVaccession[which(grepl(mim_no,clinvar_txt$PhenotypeIDS))]
  causal_variants <- unique(unlist(lapply(clinvar_accessions,retrieve_variant_consequence)))
  return(causal_variants)
}
  
universe_df$omim_causal_variant<-rep(NA,length(universe_df$omim_phen_id))

index<-0
#3185
universe_df$omim_causal_variant[which(universe_df$omim=="Y")] <-lapply(1:length(universe_df$omim_phen_id[which(universe_df$omim=="Y")]),function(x) {
  index<<-x+1
  a<-lapply(universe_df$omim_phen_id[which(universe_df$omim=="Y")][[x]],find_phenotype_causal_variants)
  cat("up to",index,"\n")
  return(a)
  })

universe_df$omim_causal_variant[which(universe_df$omim=="Y")] <-lapply(10:11,function(x) {
  index<<-x+1
  a<-lapply(universe_df$omim_phen_id[which(universe_df$omim=="Y")][[x]],find_phenotype_causal_variants)
  cat("up to",index,"\n")
  return(a)
})

#files downloaded from https://github.com/macarthur-lab/clinvar/blob/master/output/b38/single/clinvar_alleles_stats.single.b38.txt
clinvar<-read.csv("Gene_lists/ClinVar/clinvar_allele_trait_pairs.single.b38.tsv",sep="\t")



