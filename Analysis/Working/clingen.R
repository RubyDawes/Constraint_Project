clingen<-read.table(file="Gene_lists/ClinGen/ClinGen_gene_curation_list.tsv",comment.char="#",sep="\t",header=TRUE)
clingen<-clingen[,c(1,5,6)]
clingen$Gene.Symbol<-checkGeneSymbols(clingen$Gene.Symbol,unmapped.as.na=TRUE,map=hgnc.table)[[3]]


aodl_severe_hi <- read.table('Gene_lists/ClinGen/aodl_severe_hi_2015_07_20.txt')[[1]]
aodl_moderate_hi <- read.table('Gene_lists/ClinGen/aodl_moderate_hi_2015_07_20.txt')[[1]]
aodl_mild_hi <- read.table('Gene_lists/ClinGen/aodl_mild_hi_2015_08_04.txt')[[1]]

clingen$severity <- rep(NA,length(clingen$Gene.Symbol))
clingen$severity[which(clingen$Haploinsufficiency.Score=="3")] <- unlist(lapply(clingen$Gene.Symbol[which(clingen$Haploinsufficiency.Score=="3")],function(x) 
                        ifelse(x%in%aodl_severe_hi,"severe",
                           ifelse(x%in%aodl_moderate_hi,"moderate",
                                  ifelse(x%in%aodl_mild_hi,"mild",NA)))))

rm(aodl_mild_hi,aodl_moderate_hi,aodl_severe_hi)

library(rvest)

retrieve_haplo_evidence<-function(gene_name){
  url<-paste("https://www.ncbi.nlm.nih.gov/projects/dbvar/clingen/clingen_gene.cgi?sym=",gene_name,"&page=print",sep="")
  query<-read_html(url)
  text<-query %>%
    html_node("#loss_ev ul:nth-of-type(2) li") %>%
    html_text()
  return(text)
}

clingen$evidence_text<-unlist(lapply(clingen$Gene.Symbol,retrieve_haplo_evidence))

save(clingen, file="output/Data/clingen_haploinsufficiency_evidence.rda", compress="bzip2")

#dominant negative

DNorGOF<- clingen$Gene.Symbol[which(grepl("dominant-negative",clingen$evidence_text,ignore.case = TRUE)|grepl("dominant negative",clingen$evidence_text,ignore.case = TRUE)|grepl("gain of function",clingen$evidence_text,ignore.case = TRUE))]
domnegun<-universe_df[which(universe_df$gene%in%domneg),]


which(grepl("gain of function",clingen$evidence_text,ignore.case = TRUE))
