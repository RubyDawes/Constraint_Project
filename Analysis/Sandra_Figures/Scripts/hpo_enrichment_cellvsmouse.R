high_phens <- c("Abnormal cellular phenotype","Abnormality of the respiratory system","Abnormality of the skeletal system",
                "Abnormality of the endocrine system","Abnormality of the immune system","Abnormal test result","Growth abnormality", 
                "Abnormality of the integument","Abnormality of the cardiovascular system","Abnormality of blood and blood-forming tissues",
                "Abnormality of the eye","Abnormality of the musculature","Abnormality of the ear", "Abnormality of the nervous system",
                "Abnormality of connective tissue","Abnormality of metabolism/homeostasis", "Abnormality of limbs", 
                "Neoplasm", "Abnormality of the genitourinary system","Abnormality of the thoracic cavity",
                "Abnormality of the breast","Constitutional symptom","Abnormality of the voice","Abnormality of the digestive system",
                "Abnormality of head or neck","Abnormality of prenatal development or birth" )

HPO<-read.table("Gene_lists/HPO/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt",header=FALSE,sep="\t",quote="",fill=TRUE,row.names=NULL)
colnames(HPO) <- c("gene_id","gene","HPO_name","HPO_id")
hpo_ids <- unlist(lapply(high_phens,function(x) vlookup(x,HPO,lookup_column = "HPO_name",result_column = "HPO_id")))
hpo_ids <- as.character(hpo_ids)
hpo_ids[c(17,20,22,24,25)]<-c('HP:0040064','HP:0045027','HP:0025142','HP:0025031','HP:0000152')
rm(HPO)

hpo_shady <- data.frame(hpo_ids = hpo_ids, hpo_names = high_phens)
rm(hpo_ids,high_phens)

hpo_cat_count <- function(hpo_id,list_of_phens) {
  return(sum(unlist(lapply(list_of_phens, function(x) length(which(hpo_id%in%x=="TRUE"))))))
}

#combination for comparison/ checking
comb<-universe_df$hpo_ancestors[which(universe_df$omim=="Y"&((universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")|(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")))]
hpo_shady$comb_count <- unlist(lapply(hpo_shady$hpo_ids, function(x) hpo_cat_count(x,comb)))

#mouse lethal cell unessential OMIM hpo cat annotations- CATA
cata<-universe_df$hpo_ancestors[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")]
hpo_shady$cata_count <- unlist(lapply(hpo_shady$hpo_ids, function(x) hpo_cat_count(x,cata)))
hpo_shady$cata_perc <- hpo_shady$cata_count/length(cata)
  
#mouse lethal+cell essential OMIM go annotations- CATB
catb <- universe_df$hpo_ancestors[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")]
hpo_shady$catb_count <- unlist(lapply(hpo_shady$hpo_ids, function(x) hpo_cat_count(x,catb)))
hpo_shady$catb_perc <- hpo_shady$catb_count/length(catb)

hpo_shady$cat_diff_perc <- hpo_shady$cata_perc-hpo_shady$catb_perc

hpo_shady$odds_ratio <- unlist(lapply(1:length(hpo_shady$hpo_ids),function(x) {
    return(fish$estimate)
}))
hpo_shady$confint_lower <- unlist(lapply(1:length(hpo_shady$hpo_ids),function(x) {
  fish<-fisher.test(matrix(c(hpo_shady$cata_count[x],length(cata)-hpo_shady$cata_count[x],
                             hpo_shady$catb_count[x],length(catb)-hpo_shady$catb_count[x]),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
hpo_shady$confint_higher <- unlist(lapply(1:length(hpo_shady$hpo_ids),function(x) {
  fish<-fisher.test(matrix(c(hpo_shady$cata_count[x],length(cata)-hpo_shady$cata_count[x],
                             hpo_shady$catb_count[x],length(catb)-hpo_shady$catb_count[x]),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
hpo_shady$pval <- unlist(lapply(1:length(hpo_shady$hpo_ids),function(x) {
  fish<-fisher.test(matrix(c(hpo_shady$cata_count[x],length(cata)-hpo_shady$cata_count[x],
                             hpo_shady$catb_count[x],length(catb)-hpo_shady$catb_count[x]),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

#--> very little significance




