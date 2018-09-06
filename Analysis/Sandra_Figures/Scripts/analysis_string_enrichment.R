http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/114-mca-multiple-correspondence-analysis-in-r-essentials/
  
  
domen<-read.table("Gene_lists/STRING_enrichment/dominant_lethal_enrichment.csv",header=TRUE,sep=",",fill=TRUE)
recen<-read.table("Gene_lists/STRING_enrichment/recessive_lethal_enrichment.csv",header=TRUE,sep=",",fill=TRUE)


domen <- recen


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

domen <- domen[-which(domen$category=="KEGG Pathways"),]
domen$IC<-rep(NA,length(domen$term.id))
domen$IC[which(domen$category=="GO Process")] <- vlookup(domen$term.name[which(domen$category=="GO Process")],BP_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Function")] <- vlookup(domen$term.name[which(domen$category=="GO Function")],MF_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Component")] <- vlookup(domen$term.name[which(domen$category=="GO Component")],CC_IC,result_column = "IC",lookup_column = "go_term")

#domen$description[which(domen$IC<1.5)]

domen <- domen[-which(domen$IC<2|domen$IC>2.5),]

#domen <- domen[-which(domen$IC<1.5|domen$IC>3|domen$X..enriched.genes<50),]

genes<- unique(unlist(lapply(domen$enriched.genes,function(x) unlist(strsplit(as.character(x),'\\|')))))

gotab <- matrix(0, ncol = length(domen$description)+1, nrow = length(genes))
gotab <- data.frame(gotab)

names(gotab)[1]<-"genes"
gotab$genes<-genes
colnames(gotab)[2:length(gotab[1,])] <- as.character(domen$description)
#for each gene, see if it's in each GO category, 1 if yes, 0 if no
a <- lapply(domen$enriched.genes,function(x) lapply(gotab$genes,function(y) ifelse(y%in%unlist(strsplit(as.character(x),'\\|')),1,0)))
gotab[,2:length(gotab[1,])]<-lapply(2:length(gotab[1,]),function(x) unlist(a[[(x-1)]]))
gotab$genes<-factor(gotab$genes,levels=gotab$genes)
gotab[,2:length(gotab[1,])]<-lapply(2:length(gotab[1,]),function(x) ifelse(gotab[,x]==1,paste(colnames(gotab[x]),"_Y",sep=""),paste(colnames(gotab[x]),"_N",sep="")))
gotab[,2:length(gotab[1,])]<-lapply(2:length(gotab[1,]),function(x) factor(gotab[,x],levels=c(paste(colnames(gotab[x]),"_Y",sep=""),paste(colnames(gotab[x]),"_N",sep=""))))
#gotab$inheritance <- lapply(gotab$genes,function(x) ifelse((universe_df$lethal_inheritance[which(universe_df$gene==as.character(x))]=="AD"|
#                                                            universe_df$lethal_inheritance[which(universe_df$gene==as.character(x))]=="XLd"|
#                                                            universe_df$lethal_inheritance[which(universe_df$gene==as.character(x))]=="MT,XLd"),
#                                                           "Dominant","Recessive"))
#gotab[,length(gotab[1,])]<-factor(gotab[,length(gotab[1,])],levels=c("Dominant","Recessive"))
#gotab<-gotab[-which(is.na(gotab$inheritance)),]
#correspondence analysis
library(ca)
library("FactoMineR")
library("factoextra")
mca<-MCA(gotab[,2:length(gotab[1,])], ncp = 100, graph = TRUE)
fviz_mca_biplot(mca, repel = FALSE,ggtheme = theme_minimal())

fviz_mca_ind(mca, 
             habillage = "system development", # color by groups 
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, ellipse.type = "confidence",
             ggtheme = theme_minimal()) 

domen_comp <- domen[which(domen$category=="GO Component"),]
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






