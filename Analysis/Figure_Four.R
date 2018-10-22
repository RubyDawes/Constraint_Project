source("plot_functions.R")

goslim <- get_ontology("http://www.geneontology.org/ontology/go.obo", propagate_relationships = "is_a",extract_tags = "minimal")

category <- unlist(lapply(goslim$ancestors,function(x) ifelse("GO:0003674"%in%x,"MF",
                                                              ifelse("GO:0008150"%in%x,"BP",
                                                                     ifelse("GO:0005575"%in%x,"CC",NA)))))
catdf <- data.frame(id = names(category),cat = unname(category))
rm(category,goslim)

#known lethal genes####
rec<-read.table("Gene_lists/BINGO/rec_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
rec$Genes.in.test.set <- lapply(rec$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))


dom<-read.table("Gene_lists/BINGO/dom_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
dom$Genes.in.test.set <- lapply(dom$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

rec$GO.ID <- lapply(rec$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))
dom$GO.ID <- lapply(dom$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))


rec$category <- vlookup(unlist(rec$GO.ID),catdf,lookup_column = "id",result_column = "cat")
rec$category[which(is.na(rec$category))] <- c("BP","CC","CC")
dom$category <- vlookup(unlist(dom$GO.ID),catdf,lookup_column = "id",result_column = "cat")
dom$category[which(is.na(dom$category))] <- c("MF","CC","CC")




dom_lethal_CC <- dom[which(dom$category=="CC"),]
dom_lethal_CC <- dom_lethal_CC[order(dom_lethal_CC$x),]
dom_lethal_CC$inh <- "dom"
dom_lethal_CC$perc <- dom_lethal_CC$x/dom_lethal_CC$X
rec_lethal_CC <- rec[which(rec$category=="CC"),]
rec_lethal_CC <- rec_lethal_CC[order(rec_lethal_CC$x,decreasing=TRUE),]
rec_lethal_CC$inh <- "rec"
rec_lethal_CC$perc <- rec_lethal_CC$x/rec_lethal_CC$X
lethal_CC<-rbind(rec_lethal_CC,dom_lethal_CC)

vague <- c("intracellular","cell","organelle","cellular_component")
lethal_CC <- lethal_CC[-which(lethal_CC$Description%in%vague),]

lethal_CC$inh <- factor(lethal_CC$inh,levels=unique(lethal_CC$inh))

lethal_CC$no <- c(1:length(lethal_CC$GO.ID))
lethal_CC$no <- factor(lethal_CC$no,levels=lethal_CC$no)


a<-ggplot(data = lethal_CC, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(rec="#a020f0",dom="#008b45"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(lethal_CC$Description, width = 1, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))


dom_lethal_MF <- dom[which(dom$category=="MF"),]
dom_lethal_MF <- dom_lethal_MF[order(dom_lethal_MF$x),]
dom_lethal_MF$inh <- "dom"
dom_lethal_MF$perc <- dom_lethal_MF$x/dom_lethal_MF$X
rec_lethal_MF <- rec[which(rec$category=="MF"),]
rec_lethal_MF <- rec_lethal_MF[order(rec_lethal_MF$x,decreasing=TRUE),]
rec_lethal_MF$inh <- "rec"
rec_lethal_MF$perc <- rec_lethal_MF$x/rec_lethal_MF$X
lethal_MF<-rbind(rec_lethal_MF,dom_lethal_MF)

vague <- c("molecular_function","binding","protein binding","kinase activity")
lethal_MF <- lethal_MF[-which(lethal_MF$Description%in%vague),]


lethal_MF$no <- c(1:length(lethal_MF$GO.ID))
lethal_MF$no <- factor(lethal_MF$no,levels=lethal_MF$no)


b<-ggplot(data = lethal_MF, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(lethal_MF$Description, width =20, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))

dom_lethal_BP <- dom[which(dom$category=="BP"),]
dom_lethal_BP <- dom_lethal_BP[order(dom_lethal_BP$x),]
dom_lethal_BP$inh <- "dom"
dom_lethal_BP$perc <- dom_lethal_BP$x/dom_lethal_BP$X
rec_lethal_BP <- rec[which(rec$category=="BP"),]
rec_lethal_BP <- rec_lethal_BP[order(rec_lethal_BP$x,decreasing=TRUE),]
rec_lethal_BP$inh <- "rec"
rec_lethal_BP$perc <- rec_lethal_BP$x/rec_lethal_BP$X
lethal_BP<-rbind(rec_lethal_BP,dom_lethal_BP)

vague <- c("biological_process","primary metabolic process","metabolic process","regulation of biological process")
lethal_BP <- lethal_BP[-which(lethal_BP$Description%in%vague),]


lethal_BP$no <- c(1:length(lethal_BP$GO.ID))
lethal_BP$no <- factor(lethal_BP$no,levels=lethal_BP$no)

dups<-lethal_BP$Description[which(duplicated(lethal_BP$Description)=="TRUE")]
lethal_BP <- lethal_BP[-which(lethal_BP$Description%in%dups),]

c<-ggplot(data = lethal_BP, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(lethal_BP$Description, width = 25, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))

#all OMIM genes####
#write.table(universe_df$gene[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd"))],
#            "output/genelists/omim_dom.txt",row.names = FALSE,col.names=FALSE, quote = FALSE,sep="\t")

#write.table(universe_df$gene[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|
#                                                            universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR"))],
#            "output/genelists/omim_rec.txt",row.names = FALSE,col.names=FALSE, quote = FALSE,sep="\t")


recomim<-read.table("Gene_lists/BINGO/omim_rec_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
recomim$Genes.in.test.set <- lapply(recomim$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

domomim<-read.table("Gene_lists/BINGO/dom_omim_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
domomim$Genes.in.test.set <- lapply(domomim$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

recomim$GO.ID <- lapply(recomim$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))
domomim$GO.ID <- lapply(domomim$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))


recomim$category <- vlookup(unlist(recomim$GO.ID),catdf,lookup_column = "id",result_column = "cat")
recomim$category[which(is.na(recomim$category))] <- c("BP","CC","CC","CC","BP")
domomim$category <- vlookup(unlist(domomim$GO.ID),catdf,lookup_column = "id",result_column = "cat")
domomim$category[which(is.na(domomim$category))] <- c("CC","MF","BP","CC","CC")


dom_omim_CC <- domomim[which(domomim$category=="CC"),]
dom_omim_CC <- dom_omim_CC[order(dom_omim_CC$x),]
dom_omim_CC$inh <- "dom"
dom_omim_CC$perc <- dom_omim_CC$x/dom_omim_CC$X

rec_omim_CC <- recomim[which(recomim$category=="CC"),]
rec_omim_CC <- rec_omim_CC[order(rec_omim_CC$x,decreasing=TRUE),]
rec_omim_CC$inh <- "rec"
rec_omim_CC$perc <- rec_omim_CC$x/rec_omim_CC$X
omim_CC<-rbind(rec_omim_CC,dom_omim_CC)

vague <- c("intracellular","cell","organelle","cellular_component")
omim_CC <- omim_CC[-which(omim_CC$Description%in%vague),]

# taking 10 most significantly enriched terms
omim_CC <- omim_CC[-order(omim_CC$corr.p.value)[15:length(omim_CC$corr.p.value)],]

omim_CC$no <- c(1:length(omim_CC$GO.ID))
omim_CC$no <- factor(omim_CC$no,levels=omim_CC$no)

d<-ggplot(data = omim_CC, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(omim_CC$Description, width = 1, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))


dom_omim_MF <- domomim[which(domomim$category=="MF"),]
dom_omim_MF <- dom_omim_MF[order(dom_omim_MF$x),]
dom_omim_MF$inh <- "dom"
dom_omim_MF$perc <- dom_omim_MF$x/dom_omim_MF$X

rec_omim_MF <- recomim[which(recomim$category=="MF"),]
rec_omim_MF <- rec_omim_MF[order(rec_omim_MF$x,decreasing=TRUE),]
rec_omim_MF$inh <- "rec"
rec_omim_MF$perc <- rec_omim_MF$x/rec_omim_MF$X
omim_MF<-rbind(rec_omim_MF,dom_omim_MF)

vague <- c("molecular_function","binding","protein binding","kinase activity")
omim_MF <- omim_MF[-which(omim_MF$Description%in%vague),]

# taking 10 most significantly enriched terms
omim_MF <- omim_MF[-order(omim_MF$corr.p.value)[15:length(omim_MF$corr.p.value)],]

omim_MF$no <- c(1:length(omim_MF$GO.ID))
omim_MF$no <- factor(omim_MF$no,levels=omim_MF$no)

e<-ggplot(data = omim_MF, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(omim_MF$Description, width = 20, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))



dom_omim_BP <- domomim[which(domomim$category=="BP"),]
dom_omim_BP <- dom_omim_BP[order(dom_omim_BP$x),]
dom_omim_BP$inh <- "dom"
dom_omim_BP$perc <- dom_omim_BP$x/dom_omim_BP$X

rec_omim_BP <- recomim[which(recomim$category=="BP"),]
rec_omim_BP <- rec_omim_BP[order(rec_omim_BP$x,decreasing=TRUE),]
rec_omim_BP$inh <- "rec"
rec_omim_BP$perc <- rec_omim_BP$x/rec_omim_BP$X
omim_BP<-rbind(rec_omim_BP,dom_omim_BP)

vague <- c("biological_process","primary metabolic process","metabolic process","regulation of biological process")
omim_BP <- omim_BP[-which(omim_BP$Description%in%vague),]

dups<-omim_BP$Description[which(duplicated(omim_BP$Description)=="TRUE")]
omim_BP <- omim_BP[-which(omim_BP$Description%in%dups),]

# taking 10 most significantly enriched terms
omim_BP <- omim_BP[-order(omim_BP$corr.p.value)[15:length(omim_BP$corr.p.value)],]

omim_BP$no <- c(1:length(omim_BP$GO.ID))
omim_BP$no <- factor(omim_BP$no,levels=omim_BP$no)

f<-ggplot(data = omim_BP, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(omim_BP$Description, width = 25, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))



ggarrange(a,b,c,d,e,f,ncol=3,nrow=2,widths=c(1,1.2,1.3))
ggsave("output/Figures/4.pdf",height=33 , width=16, units='cm')
rm(list=setdiff(ls(),"universe_df"))
