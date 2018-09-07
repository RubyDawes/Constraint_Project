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







###bar plot of enriched BP terms in lethal dominant genes####
#####filtering enrichment terms to get good information content, significance#####
domen_process<-domen[which(domen$category=="GO Process"),]
domen_process <- domen_process[-which(domen_process$IC<3|domen_process$IC>5),]
domen_process<- domen_process[order(domen_process$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significant term####
sim<-termSim(as.character(domen_process$term.name),as.character(domen_process$term.name),bp,method = "Wang")
remove <- numeric()
for (i in 1:length(domen_process$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(domen_process$term.name%in%redundant))
}
domen_process <- domen_process[-unique(remove),]

###removing terms that have over 80% of genes in common- keeping most significant term ####

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
plot$term <- stringr::str_wrap(plot$term,25)
plot$term<-factor(plot$term,levels=plot$term)

d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#e32026",width=0.8)+bar_theme_go()
d<- d+labs(y = "Number of genes",x="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,10,20,30,40,50,60),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_bp_dominant_lethal.pdf",height=20, width=30, units='cm')








































###bar plot of enriched MF terms in lethal dominant genes####
#####filtering enrichment terms to get good information content, significance#####
domen_function<-domen[which(domen$category=="GO Function"),]
domen_function <- domen_function[-which(domen_function$IC<3.9|domen_function$IC>4.5),]
domen_function<- domen_function[order(domen_function$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significant term####
sim<-termSim(as.character(domen_function$term.name),as.character(domen_function$term.name),bp,method = "Wang")
remove <- numeric()
for (i in 1:length(domen_function$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(domen_function$term.name%in%redundant))
}
domen_function <- domen_function[-unique(remove),]


###removing terms that have over 80% of genes in common- keeping most significant termm ####

enrichgen<-lapply(domen_function$enriched.genes,function(x) unique(unlist(lapply(x,function(y)unlist(strsplit(as.character(x),'\\|')) ))))
perccommont<-lapply(1:length(enrichgen),function(y) lapply(1:length(enrichgen),function(x) length(which(enrichgen[[x]]%in%enrichgen[[y]]))/max(c(length(enrichgen[[x]]),length(enrichgen[[y]])))))
perccommon <- matrix(0, ncol = length(domen_function$description), nrow = length(domen_function$description))
perccommon <- data.frame(perccommon)
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
plot$term <- stringr::str_wrap(plot$term,25) 
plot$term<-factor(plot$term,levels=plot$term)

d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#e32026",width=0.8)+bar_theme_go()
d<- d+labs(y = "Number of genes",x="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,10,20,30,40,50,60),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_mf_dominant_lethal.pdf",height=20, width=30, units='cm')











###bar plot of enriched CC terms in lethal dominant genes####
#####filtering enrichment terms to get good information content, significance#####
domen_component<-domen[which(domen$category=="GO Component"),]
domen_component <- domen_component[-which(domen_component$IC<1.5|domen_component$IC>6),]
domen_component<- domen_component[order(domen_component$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significant term####
sim<-termSim(as.character(domen_component$term.name),as.character(domen_component$term.name),bp,method = "Wang")
remove <- numeric()
for (i in 1:length(domen_component$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(domen_component$term.name%in%redundant))
}
if (length(unique(remove))>0){
  domen_component <- domen_component[-unique(remove),]
}


###removing terms that have over 80% of genes in common- keeping most enriched term ####
enrichgen<-lapply(domen_component$enriched.genes,function(x) unique(unlist(lapply(x,function(y)unlist(strsplit(as.character(x),'\\|')) ))))
perccommont<-lapply(1:length(enrichgen),function(y) lapply(1:length(enrichgen),function(x) length(which(enrichgen[[x]]%in%enrichgen[[y]]))/max(c(length(enrichgen[[x]]),length(enrichgen[[y]])))))
perccommon <- matrix(0, ncol = length(domen_component$description), nrow = length(domen_component$description))
perccommon <- data.frame(perccommon)
rownames(perccommon)<-as.character(domen_component$term.name)
colnames(perccommon)<-as.character(domen_component$term.name)
perccommon[,1:length(colnames(perccommon))]<-lapply(1:length(colnames(perccommon)),function(x) unlist(perccommont[[x]]))

remove <- numeric()
for (i in 1:length(domen_component$category)){
  a<-which(perccommon[,i]>0.8)
  overlap<-rownames(perccommon)[a[which(a>i)]]
  remove <- append(remove,which(domen_component$term.name%in%overlap))
}
if (length(unique(remove))>0){
  domen_component <- domen_component[-unique(remove),]
}




###getting gene list from enriched terms####
domen_component<- domen_component[order(domen_component$X..enriched.genes,decreasing=TRUE),]
genes<- unique(unlist(lapply(domen_component$enriched.genes[1:7],function(x) unlist(strsplit(as.character(x),'\\|')))))

####bar plot of number of genes in enriched categories####
plot <- data.frame(term=domen_component$description[1:7])
plot$no <- vlookup(plot$term,domen_component,result_column = "X..enriched.genes",lookup_column="description")
plot$term <- stringr::str_wrap(plot$term,25) 
plot$term<-factor(plot$term,levels=plot$term)

d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#e32026",width=0.8)+bar_theme_go()
d<- d+labs(y = "Number of genes",x="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_cc_dominant_lethal.pdf",height=20, width=30, units='cm')













######loading in recessive enrichment data####
recen<-read.table("Gene_lists/STRING_enrichment/recessive_lethal_enrichment.csv",header=TRUE,sep=",",fill=TRUE)
recen <- recen[-which(recen$category=="KEGG Pathways"),]
###adding information content of enriched GO terms####
recen$IC<-rep(NA,length(recen$term.id))
recen$IC[which(recen$category=="GO Process")] <- vlookup(recen$term.name[which(recen$category=="GO Process")],BP_IC,result_column = "IC",lookup_column = "go_term")
recen$IC[which(recen$category=="GO Function")] <- vlookup(recen$term.name[which(recen$category=="GO Function")],MF_IC,result_column = "IC",lookup_column = "go_term")
recen$IC[which(recen$category=="GO Component")] <- vlookup(recen$term.name[which(recen$category=="GO Component")],CC_IC,result_column = "IC",lookup_column = "go_term")










###bar plot of enriched BP terms in lethal recessive genes####
#####filtering enrichment terms to get good information content, significance#####
recen_process<-recen[which(recen$category=="GO Process"),]
recen_process <- recen_process[-which(recen_process$IC<3|recen_process$IC>5),]
recen_process<- recen_process[order(recen_process$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significant term####
sim<-termSim(as.character(recen_process$term.name),as.character(recen_process$term.name),bp,method = "Wang")
remove <- numeric()
for (i in 1:length(recen_process$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(recen_process$term.name%in%redundant))
}
recen_process <- recen_process[-unique(remove),]

###removing terms that have over 80% of genes in common- keeping most significant term ####
enrichgen<-lapply(recen_process$enriched.genes,function(x) unique(unlist(lapply(x,function(y)unlist(strsplit(as.character(x),'\\|')) ))))
perccommont<-lapply(1:length(enrichgen),function(y) lapply(1:length(enrichgen),function(x) length(which(enrichgen[[x]]%in%enrichgen[[y]]))/max(c(length(enrichgen[[x]]),length(enrichgen[[y]])))))
perccommon <- matrix(0, ncol = length(recen_process$description), nrow = length(recen_process$description))
perccommon <- data.frame(perccommon)
rownames(perccommon)<-as.character(recen_process$term.name)
colnames(perccommon)<-as.character(recen_process$term.name)
perccommon[,1:length(colnames(perccommon))]<-lapply(1:length(colnames(perccommon)),function(x) unlist(perccommont[[x]]))

remove <- numeric()
for (i in 1:length(recen_process$category)){
  a<-which(perccommon[,i]>0.8)
  overlap<-rownames(perccommon)[a[which(a>i)]]
  remove <- append(remove,which(recen_process$term.name%in%overlap))
}
recen_process <- recen_process[-unique(remove),]




###getting gene list from enriched terms####
recen_process<- recen_process[order(recen_process$X..enriched.genes,decreasing=TRUE),]
genes<- unique(unlist(lapply(recen_process$enriched.genes[1:7],function(x) unlist(strsplit(as.character(x),'\\|')))))

####bar plot of number of genes in enriched categories####
plot <- data.frame(term=recen_process$description[1:7])
plot$no <- vlookup(plot$term,recen_process,result_column = "X..enriched.genes",lookup_column="description")
plot$term <- stringr::str_wrap(plot$term,25) 
plot$term<-factor(plot$term,levels=plot$term)

d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#2c5daa",width=0.8)+bar_theme_go()
d<- d+labs(y = "Number of genes",x="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,30,60,90,120,150),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_bp_recessive_lethal.pdf",height=20, width=30, units='cm')




###bar plot of enriched MF terms in lethal dominant genes####
#####filtering enrichment terms to get good information content, significance#####
recen_function<-recen[which(recen$category=="GO Function"),]
recen_function <- recen_function[-which(recen_function$IC<2.5|recen_function$IC>5),]
recen_function<- recen_function[order(recen_function$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significant term####
sim<-termSim(as.character(recen_function$term.name),as.character(recen_function$term.name),bp,method = "Wang")
remove <- numeric()
for (i in 1:length(recen_function$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(recen_function$term.name%in%redundant))
}
recen_function <- recen_function[-unique(remove),]


###removing terms that have over 80% of genes in common- keeping most significant termm ####

enrichgen<-lapply(recen_function$enriched.genes,function(x) unique(unlist(lapply(x,function(y)unlist(strsplit(as.character(x),'\\|')) ))))
perccommont<-lapply(1:length(enrichgen),function(y) lapply(1:length(enrichgen),function(x) length(which(enrichgen[[x]]%in%enrichgen[[y]]))/max(c(length(enrichgen[[x]]),length(enrichgen[[y]])))))
perccommon <- matrix(0, ncol = length(recen_function$description), nrow = length(recen_function$description))
perccommon <- data.frame(perccommon)
rownames(perccommon)<-as.character(recen_function$term.name)
colnames(perccommon)<-as.character(recen_function$term.name)
perccommon[,1:length(colnames(perccommon))]<-lapply(1:length(colnames(perccommon)),function(x) unlist(perccommont[[x]]))

remove <- numeric()
for (i in 1:length(recen_function$category)){
  a<-which(perccommon[,i]>0.8)
  overlap<-rownames(perccommon)[a[which(a>i)]]
  remove <- append(remove,which(recen_function$term.name%in%overlap))
}
recen_function <- recen_function[-unique(remove),]




###getting gene list from enriched terms####
recen_function<- recen_function[order(recen_function$X..enriched.genes,decreasing=TRUE),]
genes<- unique(unlist(lapply(recen_function$enriched.genes[1:7],function(x) unlist(strsplit(as.character(x),'\\|')))))

####bar plot of number of genes in enriched categories####
plot <- data.frame(term=recen_function$description[1:7])
plot$no <- vlookup(plot$term,recen_function,result_column = "X..enriched.genes",lookup_column="description")
plot$term <- stringr::str_wrap(plot$term,25) 
plot$term<-factor(plot$term,levels=plot$term)

d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#2c5daa",width=0.8)+bar_theme_go()
d<- d+labs(y = "Number of genes",x="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,30,60,90,120,150),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_mf_recessive_lethal.pdf",height=20, width=30, units='cm')














###bar plot of enriched CC terms in lethal dominant genes####
#####filtering enrichment terms to get good information content, significance#####
recen_component<-recen[which(recen$category=="GO Component"),]
recen_component <- recen_component[-which(recen_component$IC<1.5|recen_component$IC>6),]
recen_component<- recen_component[order(recen_component$FDR.p.value),]

####semantic similarity of enriched GO terms- removing redundant terms, keeping most significant term####
sim<-termSim(as.character(recen_component$term.name),as.character(recen_component$term.name),bp,method = "Wang")
remove <- numeric()
for (i in 1:length(recen_component$category)){
  a<-which(sim[,i]>0.7)
  redundant<-names(a[which(a>i)])
  remove <- append(remove,which(recen_component$term.name%in%redundant))
}
if (length(unique(remove))>0){
  recen_component <- recen_component[-unique(remove),]
}


###removing terms that have over 80% of genes in common- keeping most enriched term ####
enrichgen<-lapply(recen_component$enriched.genes,function(x) unique(unlist(lapply(x,function(y)unlist(strsplit(as.character(x),'\\|')) ))))
perccommont<-lapply(1:length(enrichgen),function(y) lapply(1:length(enrichgen),function(x) length(which(enrichgen[[x]]%in%enrichgen[[y]]))/max(c(length(enrichgen[[x]]),length(enrichgen[[y]])))))
perccommon <- matrix(0, ncol = length(recen_component$description), nrow = length(recen_component$description))
perccommon <- data.frame(perccommon)
rownames(perccommon)<-as.character(recen_component$term.name)
colnames(perccommon)<-as.character(recen_component$term.name)
perccommon[,1:length(colnames(perccommon))]<-lapply(1:length(colnames(perccommon)),function(x) unlist(perccommont[[x]]))

remove <- numeric()
for (i in 1:length(recen_component$category)){
  a<-which(perccommon[,i]>0.8)
  overlap<-rownames(perccommon)[a[which(a>i)]]
  remove <- append(remove,which(recen_component$term.name%in%overlap))
}
if (length(unique(remove))>0){
  recen_component <- recen_component[-unique(remove),]
}




###getting gene list from enriched terms####
recen_component<- recen_component[order(recen_component$X..enriched.genes,decreasing=TRUE),]
genes<- unique(unlist(lapply(recen_component$enriched.genes[1:7],function(x) unlist(strsplit(as.character(x),'\\|')))))

####bar plot of number of genes in enriched categories####
plot <- data.frame(term=recen_component$description[1:7])
plot$no <- vlookup(plot$term,recen_component,result_column = "X..enriched.genes",lookup_column="description")
plot$term <- stringr::str_wrap(plot$term,25)
plot$term<-factor(plot$term,levels=plot$term)
d <- ggplot(dat=plot, aes(x=term,y=no))
d<- d+geom_bar(stat="identity",fill = "#2c5daa",width=0.8)+bar_theme_go()
d<- d+labs(y = "Number of genes",x="Enriched GO BP term")
d<- d+scale_y_continuous(breaks = c(0,20,40,60,80,100,120,140,160),expand = c(0, 0))+coord_flip()
d
ggsave("Analysis/Sandra_Figures/Figs/gobar_cc_recessive_lethal.pdf",height=20, width=30, units='cm')







