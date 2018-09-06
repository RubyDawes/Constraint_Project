#all dominant lethal genes
domgen<-universe_df$gene[which(universe_df$human_lethal=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd"))]
#all recessive lethal genes
recgen<-universe_df$gene[which(universe_df$human_lethal=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="XLr"))]


###information content from GOSemSim######
library(GOSemSim)
#source("https://bioconductor.org/biocLite.R")
#biocLite("org.Hs.eg.db")
OrgDb<-load_OrgDb("org.Hs.eg.db")
bp<-godata(OrgDb = OrgDb, keytype = "ENTREZID", 'BP', computeIC = TRUE)
BP_IC <- data.frame(go_term= names(bp@IC), IC = as.vector(bp@IC))
mf<-godata(OrgDb = OrgDb, keytype = "ENTREZID", 'MF', computeIC = TRUE)
MF_IC <- data.frame(go_term= names(mf@IC), IC = as.vector(mf@IC))
cc<-godata(OrgDb = OrgDb, keytype = "ENTREZID", 'CC', computeIC = TRUE)
CC_IC <- data.frame(go_term= names(cc@IC), IC = as.vector(cc@IC))

######loading in enrichment data####
domen<-read.table("Gene_lists/STRING_enrichment/dominant_lethal_enrichment.csv",header=TRUE,sep=",",fill=TRUE)
domen <- domen[-which(domen$category=="KEGG Pathways"),]

###adding information content of enriched GO terms####
domen$IC<-rep(NA,length(domen$term.id))
domen$IC[which(domen$category=="GO Process")] <- vlookup(domen$term.name[which(domen$category=="GO Process")],BP_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Function")] <- vlookup(domen$term.name[which(domen$category=="GO Function")],MF_IC,result_column = "IC",lookup_column = "go_term")
domen$IC[which(domen$category=="GO Component")] <- vlookup(domen$term.name[which(domen$category=="GO Component")],CC_IC,result_column = "IC",lookup_column = "go_term")



#bar plot of enriched BP terms in lethal dominant genes
#####filtering enrichment terms to get good information content, significance#####
domen_process<-domen[which(domen$category=="GO Process"),]
domen_process <- domen_process[-which(domen_process$IC<3|domen_process$IC>5),]
domen_process<- domen_process[order(domen_process$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significantly enriched term####
library(GOSemSim)
sim<-termSim(as.character(domen_process$term.name),as.character(domen_process$term.name),bp,method = "Wang")

remove <- numeric()
for (i in 1:length(domen_process$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(domen_process$term.name%in%redundant))
}
domen_process <- domen_process[-unique(remove),]


###removing terms that have over 80% of genes in common- keeping most significantly enriched term ####
enrichgen<-lapply(domen_process$enriched.genes,function(x) unique(unlist(lapply(x,function(y)unlist(strsplit(as.character(x),'\\|')) ))))
perccommont<-lapply(1:length(enrichgen),function(y) lapply(1:length(enrichgen),function(x) length(which(enrichgen[[x]]%in%enrichgen[[y]]))/max(c(length(enrichgen[[x]]),length(enrichgen[[y]])))))
perccommon <- matrix(0, ncol = length(domen_process$description), nrow = length(domen_process$description))
perccommon <- data.frame(perccommon)
rownames(perccommon)<-as.character(domen_process$term.name)
colnames(perccommon)<-as.character(domen_process$term.name)
perccommon[,1:length(colnames(perccommon))]<-lapply(1:length(colnames(perccommon)),function(x) unlist(perccommont[[x]]))

remove <- numeric()
for (i in 1:length(domen_process$category)){
  a<-which(perccommon[,i]>0.8)
  overlap<-rownames(perccommon)[a[which(a>i)]]
  remove <- append(remove,which(domen_process$term.name%in%overlap))
}
domen_process <- domen_process[-unique(remove),]




###getting gene list from enriched terms####
domen_process<- domen_process[order(domen_process$X..enriched.genes,decreasing=TRUE),]
genes<- unique(unlist(lapply(domen_process$enriched.genes[1:7],function(x) unlist(strsplit(as.character(x),'\\|')))))

####bar plot of number of genes in enriched categories####
plot <- data.frame(term=domen_process$description[1:7])
plot$no <- vlookup(plot$term,domen_process,result_column = "X..enriched.genes",lookup_column="description")
plot$term<-factor(plot$term,levels=plot$term)

d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#e32026",width=0.8)+bar_theme_go()
d<- d+labs(x = "Number of genes",y="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,10,20,30,40,50,60),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_bp_dominant_lethal.pdf",height=20, width=20, units='cm')

























####looking for 5 enriched terms which maximise number of represented genes####
y<-1000
x=data.frame(numbers=I(as.list(rep(NA,y))),gene_no=rep(NA,y))

for (i in 1:y) {
  num <- sample(1:length(domen_process$term.id),10,replace=F)
  x$numbers[[i]]<-num
  x$gene_no[[i]]<-length(unique(unlist(lapply(domen_process$enriched.genes[num],function(x) unlist(strsplit(as.character(x),'\\|'))))))
}

max(x$gene_no)
x$numbers[which(x$gene_no==max(x$gene_no))]
domen_process$description[unlist(x$numbers[which(x$gene_no==max(x$gene_no))][[1]])]















#bar plot of enriched MF terms in lethal dominant genes
#####filtering enrichment terms to get good information content, significance#####
domen_function<-domen[which(domen$category=="GO Function"),]
domen_function <- domen_function[-which(domen_function$IC<2|domen_function$IC>4.5),]
domen_function<- domen_function[order(domen_function$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significantly enriched term####

sim<-termSim(as.character(domen_function$term.name),as.character(domen_function$term.name),bp,method = "Wang")

remove <- numeric()
for (i in 1:length(domen_function$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(domen_function$term.name%in%redundant))
}
domen_function <- domen_function[-unique(remove),]


###removing terms that have over 80% of genes in common- keeping most significantly enriched term ####
enrichgen<-lapply(domen_function$enriched.genes,function(x) unique(unlist(lapply(x,function(y)unlist(strsplit(as.character(x),'\\|')) ))))
perccommont<-lapply(1:length(enrichgen),function(y) lapply(1:length(enrichgen),function(x) length(which(enrichgen[[x]]%in%enrichgen[[y]]))/max(c(length(enrichgen[[x]]),length(enrichgen[[y]])))))
perccommon <- matrix(0, ncol = length(domen_function$description), nrow = length(domen_function$description))
perccommon <- data.frame(gotab)
rownames(perccommon)<-as.character(domen_function$term.name)
colnames(perccommon)<-as.character(domen_function$term.name)
perccommon[,1:length(colnames(perccommon))]<-lapply(1:length(colnames(perccommon)),function(x) unlist(perccommont[[x]]))

remove <- numeric()
for (i in 1:length(domen_function$category)){
  a<-which(perccommon[,i]>0.8)
  overlap<-rownames(perccommon)[a[which(a>i)]]
  remove <- append(remove,which(domen_function$term.name%in%overlap))
}
domen_function <- domen_function[-unique(remove),]




###getting gene list from enriched terms####
domen_function<- domen_function[order(domen_function$X..enriched.genes,decreasing=TRUE),]
genes<- unique(unlist(lapply(domen_function$enriched.genes[1:7],function(x) unlist(strsplit(as.character(x),'\\|')))))

####bar plot of number of genes in enriched categories####
plot <- data.frame(term=domen_function$description[1:7])
plot$no <- vlookup(plot$term,domen_function,result_column = "X..enriched.genes",lookup_column="description")
plot$term<-factor(plot$term,levels=plot$term)

d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#e32026",width=0.8)+bar_theme_go()
d<- d+labs(x = "Number of genes",y="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,10,20,30,40,50,60),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_bp_dominant_lethal.pdf",height=20, width=20, units='cm')






####looking for 5 enriched terms which maximise number of represented genes####
y<-1000
x=data.frame(numbers=I(as.list(rep(NA,y))),gene_no=rep(NA,y))

for (i in 1:y) {
  num <- sample(1:length(domen_function$term.id),10,replace=F)
  x$numbers[[i]]<-num
  x$gene_no[[i]]<-length(unique(unlist(lapply(domen_function$enriched.genes[num],function(x) unlist(strsplit(as.character(x),'\\|'))))))
}

max(x$gene_no)
x$numbers[which(x$gene_no==max(x$gene_no))]
domen_function$description[unlist(x$numbers[which(x$gene_no==max(x$gene_no))][[1]])]

















