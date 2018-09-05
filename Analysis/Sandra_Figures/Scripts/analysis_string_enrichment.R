domen<-read.table("Gene_lists/STRING_enrichment/dominant_lethal_enrichment.csv",header=TRUE,sep=",",fill=TRUE)
recen<-read.table("Gene_lists/STRING_enrichment/recessive_lethal_enrichment.csv",header=TRUE,sep=",",fill=TRUE)

library(GOSemSim)
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
OrgDb<-load_OrgDb("org.Hs.eg.db")
a<-godata(OrgDb = OrgDb, keytype = "ENTREZID", 'BP', computeIC = TRUE)
BP_IC <- data.frame(go_term= names(a@IC), IC = as.vector(a@IC))
a<-godata(OrgDb = OrgDb, keytype = "ENTREZID", 'MF', computeIC = TRUE)
MF_IC <- data.frame(go_term= names(a@IC), IC = as.vector(a@IC))
a<-godata(OrgDb = OrgDb, keytype = "ENTREZID", 'CC', computeIC = TRUE)
CC_IC <- data.frame(go_term= names(a@IC), IC = as.vector(a@IC))

domen$IC<-rep(NA,length(domen$term.id))
domen$IC[which(domen$category=="GO Process")] <- vlookup(domen$term.name[which(domen$category=="GO Process")],BP_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Function")] <- vlookup(domen$term.name[which(domen$category=="GO Function")],MF_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Component")] <- vlookup(domen$term.name[which(domen$category=="GO Component")],CC_IC,result_column = "IC",lookup_column = "go_term")




domen$description[which(domen$category=="GO Process")][order(domen$IC[which(domen$category=="GO Process")])[1:5]]
domen$description[which(domen$category=="GO Process")][order(domen$IC[which(domen$category=="GO Process")])[1000:1120]]

domen$description[which(domen$category=="GO Process"&domen$IC>8&domen$IC<10)]
domen$description[which(domen$category=="GO Process"&domen$IC>=0&domen$IC<1)]

domen$description[which(domen$category=="GO Process"&domen$IC>1&domen$IC<2)]

domen$description[which(domen$category=="GO Process"&domen$IC>2&domen$IC<3&domen$FDR.p.value<0.005)]
domen$description[which(domen$category=="GO Process"&domen$IC>3&domen$IC<4)]
domen$description[which(domen$category=="GO Process"&domen$IC>4&domen$IC<6)]
domen$description[which(domen$category=="GO Process"&domen$IC>6&domen$IC<8)]



domen$description[which(domen$category=="GO Process"&domen$IC>2&domen$IC<3&domen$FDR.p.value<0.005&domen$X..enriched.genes>70)]






