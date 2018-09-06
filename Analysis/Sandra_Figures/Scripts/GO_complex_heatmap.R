
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

domen<-read.table("Gene_lists/STRING_enrichment/dominant_lethal_enrichment.csv",header=TRUE,sep=",",fill=TRUE)

domen <- domen[-which(domen$category=="KEGG Pathways"),]
domen$IC<-rep(NA,length(domen$term.id))
domen$IC[which(domen$category=="GO Process")] <- vlookup(domen$term.name[which(domen$category=="GO Process")],BP_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Function")] <- vlookup(domen$term.name[which(domen$category=="GO Function")],MF_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Component")] <- vlookup(domen$term.name[which(domen$category=="GO Component")],CC_IC,result_column = "IC",lookup_column = "go_term")
domen <- domen[-which(domen$IC<1.5|domen$IC>3),]

  domen <- domen[which(domen$category=="GO Component"),]

genes<- unique(unlist(lapply(domen$enriched.genes,function(x) unlist(strsplit(as.character(x),'\\|')))))

gotab <- matrix(0, ncol = length(domen$description)+1, nrow = length(genes))
gotab <- data.frame(gotab)


names(gotab)[1]<-"genes"
gotab$genes<-genes
colnames(gotab)[2:length(gotab[1,])] <- as.character(domen$description)
#for each gene, see if it's in each GO category, 1 if yes, 0 if no
a <- lapply(domen$enriched.genes,function(x) lapply(gotab$genes,function(y) ifelse(y%in%unlist(strsplit(as.character(x),'\\|')),1,0)))
gotab[,2:length(gotab[1,])]<-lapply(2:length(gotab[1,]),function(x) unlist(a[[(x-1)]]))

gotab.m <- melt(gotab)

ggplot(data = gotab.m, aes(x = genes, y = variable)) +
  geom_tile(aes(fill = value))+ scale_fill_gradient(low = "white",high = "steelblue")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

rownames(gotab)<-gotab[,1]
png(file = "heatmap2.png",width = 50, height = 50, units = "cm",res=1000)
heatmap.2(data.matrix(gotab[,-1]),dendrogram="none",keysize=0,trace="none")
dev.off()
ggsave("heatmap.pdf",height=50, width=50, units='cm')

