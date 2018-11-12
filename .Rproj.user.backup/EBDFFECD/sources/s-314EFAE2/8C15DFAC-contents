
library(STRINGdb)

lethal <- universe_df[which(universe_df$human_lethal_B=="Y"),c("gene","lethal_inheritance")]

string_db <- STRINGdb$new( version="10", species=9606, 
              score_threshold=0, input_directory="" )
background <-string_db$map(universe_df,"gene",removeUnmappedRows = TRUE)
string_db$set_background(background$STRING_id)
save(string_db,file="output/Data/string_db.rda",compress="bzip2")

a<-string_db$map(nucleus,"V1")
test<-string_db$get_enrichment(a$STRING_id,category = "Component", methodMT = "bonferroni", iea = TRUE )
test<-string_db$get_enrichment(background$STRING_id, category = "Process", methodMT = "bonferroni", iea = TRUE )
testCC<-string_db$get_enrichment(background$STRING_id, category = "Component", methodMT = "bonferroni", iea = TRUE )

lethal_mapped <- string_db$map( lethal, "gene", removeUnmappedRows = FALSE )
lethal_mapped$domrec <- ifelse((lethal_mapped$lethal_inheritance=="AR"|lethal_mapped$lethal_inheritance=="XLr"|lethal_mapped$lethal_inheritance=="MT,AR"),"recessive",
                               ifelse((lethal_mapped$lethal_inheritance=="AD"|lethal_mapped$lethal_inheritance=="XLd"|lethal_mapped$lethal_inheritance=="MT,XLd"),"dominant",NA))


rec_enr_pr<-string_db$get_enrichment( lethal_mapped[which(lethal_mapped$domrec=="recessive"&!is.na(lethal_mapped$STRING_id)),"STRING_id"], 
                                   category = "Process", methodMT = "bonferroni", iea = TRUE )
rec_enr_pr$category <- rep("Process",length(rec_enr_pr$term_id))
rec_enr_fn<-string_db$get_enrichment( lethal_mapped[which(lethal_mapped$domrec=="recessive"&!is.na(lethal_mapped$STRING_id)),"STRING_id"], 
                                      category = "Function", methodMT = "bonferroni", iea = TRUE )
rec_enr_fn$category <- rep("Function",length(rec_enr_fn$term_id))
rec_enr_co<-string_db$get_enrichment( lethal_mapped[which(lethal_mapped$domrec=="recessive"&!is.na(lethal_mapped$STRING_id)),"STRING_id"], 
                                      category = "Component", methodMT = "bonferroni", iea = TRUE )
rec_enr_co$category <- rep("Component",length(rec_enr_co$term_id))
rec_enr<-rbind(rec_enr_pr,rec_enr_fn,rec_enr_co)
#rm(rec_enr_pr,rec_enr_fn,rec_enr_co)

dom_enr_pr<-string_db$get_enrichment( lethal_mapped[which(lethal_mapped$domrec=="dominant"&!is.na(lethal_mapped$STRING_id)),"STRING_id"], 
                                      category = "Process", methodMT = "bonferroni", iea = TRUE )
dom_enr_pr$category <- rep("Process",length(dom_enr_pr$term_id))
dom_enr_fn<-string_db$get_enrichment( lethal_mapped[which(lethal_mapped$domrec=="dominant"&!is.na(lethal_mapped$STRING_id)),"STRING_id"], 
                                      category = "Function", methodMT = "bonferroni", iea = TRUE )
dom_enr_fn$category <- rep("Function",length(dom_enr_fn$term_id))
dom_enr_co<-string_db$get_enrichment( lethal_mapped[which(lethal_mapped$domrec=="dominant"&!is.na(lethal_mapped$STRING_id)),"STRING_id"], 
                                      category = "Component", methodMT = "bonferroni", iea = FALSE)
dom_enr_co$category <- rep("Component",length(dom_enr_co$term_id))
dom_enr<-rbind(dom_enr_pr,dom_enr_fn,dom_enr_co)
#rm(dom_enr_pr,dom_enr_fn,dom_enr_co)

rec_enr <- rec_enr[which(rec_enr$pvalue_fdr<0.05),]
dom_enr <- dom_enr[which(dom_enr$pvalue_fdr<0.05),]

retrieve_enriched_genes <- function(go_id,inheritance) {
  prot <- string_db$get_term_proteins(go_id,enableIEA=TRUE)
  prot <- prot[which(prot$preferred_name%in%lethal_mapped[which(lethal_mapped$domrec==inheritance&!is.na(lethal_mapped$STRING_id)),"gene"]),]
  return(prot$preferred_name)
}
rec_prot<- lapply(rec_enr$term_id,function(x) retrieve_enriched_genes(x,"recessive"))
rec_enr$genes <- rec_prot
rec_enr$lengths <- lengths(rec_enr$genes)
save(rec_enr,file="output/Data/recessive_lethal_enrichment.rda",compress="bzip2")
dom_prot<- lapply(dom_enr$term_id,function(x) retrieve_enriched_genes(x,"dominant"))
dom_enr$genes <- dom_prot
dom_enr$lengths <- lengths(dom_enr$genes)
save(dom_enr,file="output/Data/dominant_lethal_enrichment.rda",compress="bzip2")


dom_enr <- dom_enr[order(dom_enr$pvalue_bonferroni),]


library(GO.db)
offspring = c(GOMFOFFSPRING[["GO:0003676"]],GOMFOFFSPRING[["GO:0140110"]])
dna_binding_terms <- dom_enr$term_description[which(dom_enr$term_id%in%offspring)]

dom_enr$term_description[which(dom_enr$term_id%in%GOCCOFFSPRING[["GO:0005634"]])]

nucprot<-string_db$get_term_proteins("GO:0005634",enableIEA=TRUE)
nucprot <- nucprot[which(nucprot$preferred_name%in%lethal_mapped[which(lethal_mapped$domrec=="recessive"&!is.na(lethal_mapped$STRING_id)),"gene"]),]

nucprot <- nucprot[which(nucprot$preferred_name%in%lethal_mapped[which(lethal_mapped$domrec=="dominant"&!is.na(lethal_mapped$STRING_id)),"gene"]),]
length(which(lethal_mapped$domrec=="dominant"&!is.na(lethal_mapped$STRING_id)))

lethal_enr <- string_db$get_enrichment( lethal_mapped$STRING_id, 
                                        category = "Component", methodMT = "bonferroni", iea = TRUE,minScore=NULL)


library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
a<-GOCollection(unlist(universe_df$go_terms[[2]]))
goSlim(a, slim, "MF")
goSlim(a, slim, "BP")
goSlim(a, slim, "CC")


goSlim(GOCollection("GO:0005634"),slim,"CC")

lapply(universe_df$go_terms[[2]],function(x) GOMFANCESTOR[[x]])


nucprot<-string_db$get_term_proteins(c("GO:0005634",GOCCOFFSPRING[["GO:0005634"]]),enableIEA=TRUE)
nucprot <- nucprot[which(nucprot$STRING_id%in%lethal_mapped[which(lethal_mapped$domrec=="dominant"&!is.na(lethal_mapped$STRING_id)),"STRING_id"]),]

  
  