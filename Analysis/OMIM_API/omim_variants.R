

extract_phen_id <- function(phenotype){
  no_phens <- length(phenotype)
  start<- unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-7))
  end<-unlist(lapply(1:no_phens,function(x) regexpr("\\(3\\)",phenotype[x])[1]-2))
  return(substr(phenotype,start,end))
}

universe_df$omim_phen_id<-lapply(universe_df$phenotype,extract_phen_id)

clinvar_txt <- read.csv("Gene_lists/ClinVar/variant_summary.txt",sep="\t",stringsAsFactors = FALSE)
clinvar_txt<- clinvar_txt[which(clinvar_txt$ClinSigSimple==1&clinvar_txt$Assembly=="GRCh38"),]

#files downloaded from https://github.com/macarthur-lab/clinvar/blob/master/output/b38/single/clinvar_alleles_stats.single.b38.txt
clinvar_molcon<-read.csv("Gene_lists/ClinVar/clinvar_allele_trait_pairs.single.b38.tsv",sep="\t")


clinvar_molcon$molecular_consequence_summary <-lapply(1:length(clinvar_molcon$molecular_consequence),function(x) unlist(strsplit(as.character(clinvar_molcon$molecular_consequence[x]),":"))[3])
clinvar_molcon$molecular_consequence_summary <-lapply(1:length(clinvar_molcon$molecular_consequence_summary),function(x) unlist(strsplit(as.character(clinvar_molcon$molecular_consequence_summary[x]),";"))[1])


phen_mim_no <- universe_df$omim_phen_id[which(universe_df$omim=="Y")][2]

find_phenotype_causal_variants <- function(phen_mim_no){
  no_nos <- length(phen_mim_no[[1]])
  mim_no <- unlist(lapply(1:no_nos,function(x) paste("OMIM:",phen_mim_no[[1]][x],sep="")))
  clinvar_accessions <- lapply(1:no_nos,function(x) clinvar_txt$RCVaccession[which(grepl(mim_no[x],clinvar_txt$PhenotypeIDS))])
  clinvar_accessions<- lapply(1:no_nos,function(x) unlist(strsplit(paste(clinvar_accessions[[x]],collapse=";"),";")))
  causal_variants <- lapply(1:no_nos,function(x) unique(unlist(clinvar_molcon$molecular_consequence_summary[which(clinvar_molcon$rcv%in%clinvar_accessions[[x]])])))
  return(causal_variants)
}

phenotype_causal_variants<-lapply(1:length(universe_df$omim_phen_id[which(universe_df$omim=="Y")]),
                                                                            function(x){
                                                                              cat("up to", x,"\n")
                                                                              a<-find_phenotype_causal_variants(universe_df$omim_phen_id[which(universe_df$omim=="Y")][x])
                                                                              return(a)
                                                                            })
universe_df$phenotype_causal_variants<-rep(NA,length(universe_df$omim_phen_id))
universe_df$phenotype_causal_variants[which(universe_df$omim=="Y")]<-phenotype_causal_variants
  
universe_df$phenotype_causal_variants_simp<-rep(NA,length(universe_df$omim_phen_id))
universe_df$phenotype_causal_variants_simp <- lapply(1:length(universe_df$phenotype_causal_variants),function(x) unique(unlist(universe_df$phenotype_causal_variants[x])))

universe_df$phenotype[which(universe_df$phenotype_causal_variants_simp=="missense variant"&universe_df$mis_z>=3.09&universe_df$Inheritance_pattern=="AD")]
which(universe_df$phenotype_causal_variants_simp=="missense variant"&universe_df$mis_z<3.09)
which(universe_df$phenotype_causal_variants_simp=="nonsense"&universe_df$pLI>=0.9)
universe_df$hpo_names[which(universe_df$phenotype_causal_variants_simp=="nonsense"&universe_df$pLI<0.9&universe_df$Inheritance_pattern=="AD")]
universe_df$hpo_names[which(universe_df$phenotype_causal_variants_simp=="nonsense"&universe_df$pLI>=0.9&universe_df$Inheritance_pattern=="AD")]



