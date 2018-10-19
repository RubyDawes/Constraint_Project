goslim <- get_ontology("Gene_lists/GOSLIM/go.obo", propagate_relationships = "is_a",extract_tags = "minimal")

category <- unlist(lapply(goslim$ancestors,function(x) ifelse("GO:0003674"%in%x,"MF",
                                                              ifelse("GO:0008150"%in%x,"BP",
                                                                     ifelse("GO:0005575"%in%x,"CC",NA)))))
catdf <- data.frame(id = names(category),cat = unname(category))
rm(category,goslim)

#known lethal genes####
rec<-read.table("Gene_lists/BINGO/rec_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
rec$Genes.in.test.set <- lapply(rec$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))
rec$category <- 
  
dom<-read.table("Gene_lists/BINGO/dom_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
dom$Genes.in.test.set <- lapply(dom$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

rec$GO.ID <- lapply(rec$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))
dom$GO.ID <- lapply(dom$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))


rec$category <- vlookup(unlist(rec$GO.ID),catdf,lookup_column = "id",result_column = "cat")
rec$category[which(is.na(rec$category))] <- c("BP","CC","CC")
dom$category <- vlookup(unlist(dom$GO.ID),catdf,lookup_column = "id",result_column = "cat")
dom$category[which(is.na(dom$category))] <- c("MF","CC","CC")




dom_lethal_CC <- dom[which(dom$category=="CC"),]
dom_lethal_CC <- dom_lethal_CC[order(dom_lethal_CC$x,decreasing=TRUE),]
dom_lethal_CC$inh <- "dom"
dom_lethal_CC$perc <- dom_lethal_CC$x/dom_lethal_CC$X
rec_lethal_CC <- rec[which(rec$category=="CC"),]
rec_lethal_CC <- rec_lethal_CC[order(rec_lethal_CC$x),]
rec_lethal_CC$inh <- "rec"
rec_lethal_CC$perc <- rec_lethal_CC$x/rec_lethal_CC$X
lethal_CC<-rbind(dom_lethal_CC,rec_lethal_CC)

vague <- c("intracellular","cell","organelle","cellular_component")
lethal_CC <- lethal_CC[-which(lethal_CC$Description%in%vague),]


lethal_CC$no <- c(1:length(lethal_CC$GO.ID))
lethal_CC$no <- factor(lethal_CC$no,levels=lethal_CC$no)

a<-ggplot(data = lethal_CC, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
    scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(lethal_CC$Description, width = 1, simplify = FALSE), paste, collapse="\n"))+
    scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))

dom_lethal_MF <- dom[which(dom$category=="MF"),]
dom_lethal_MF <- dom_lethal_MF[order(dom_lethal_MF$x,decreasing=TRUE),]
dom_lethal_MF$inh <- "dom"
dom_lethal_MF$perc <- dom_lethal_MF$x/dom_lethal_CC$X
rec_lethal_MF <- rec[which(rec$category=="MF"),]
rec_lethal_MF <- rec_lethal_MF[order(rec_lethal_MF$x),]
rec_lethal_MF$inh <- "rec"
rec_lethal_MF$perc <- rec_lethal_MF$x/rec_lethal_MF$X
lethal_MF<-rbind(dom_lethal_MF,rec_lethal_MF)

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
dom_lethal_BP <- dom_lethal_BP[order(dom_lethal_BP$x,decreasing=TRUE),]
dom_lethal_BP$inh <- "dom"
dom_lethal_BP$perc <- dom_lethal_BP$x/dom_lethal_CC$X
rec_lethal_BP <- rec[which(rec$category=="BP"),]
rec_lethal_BP <- rec_lethal_BP[order(rec_lethal_BP$x),]
rec_lethal_BP$inh <- "rec"
rec_lethal_BP$perc <- rec_lethal_BP$x/rec_lethal_BP$X
lethal_BP<-rbind(dom_lethal_BP,rec_lethal_BP)

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
write.table(universe_df$gene[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd"))],
            "output/genelists/omim_dom.txt",row.names = FALSE,col.names=FALSE, quote = FALSE,sep="\t")

write.table(universe_df$gene[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|
                                                            universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR"))],
            "output/genelists/omim_rec.txt",row.names = FALSE,col.names=FALSE, quote = FALSE,sep="\t")


recomim<-read.table("Gene_lists/BINGO/omim_rec_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
recomim$Genes.in.test.set <- lapply(recomim$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))
  
domomim<-read.table("Gene_lists/BINGO/dom_omim_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
domomim$Genes.in.test.set <- lapply(domomim$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

recomim$GO.ID <- lapply(recomim$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))
domomim$GO.ID <- lapply(domomim$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))


recomim$category <- vlookup(unlist(recomim$GO.ID),catdf,lookup_column = "id",result_column = "cat")
recomim$category[which(is.na(recomim$category))] <- c("MF","CC","CC","CC","BP")
domomim$category <- vlookup(unlist(domomim$GO.ID),catdf,lookup_column = "id",result_column = "cat")
domomim$category[which(is.na(domomim$category))] <- c("CC","MF","BP","CC","CC")


dom_omim_CC <- domomim[which(domomim$category=="CC"),]
dom_omim_CC <- dom_omim_CC[order(dom_omim_CC$x,decreasing=TRUE),]
dom_omim_CC$inh <- "dom"
dom_omim_CC$perc <- dom_omim_CC$x/dom_omim_CC$X

rec_omim_CC <- recomim[which(recomim$category=="CC"),]
rec_omim_CC <- rec_omim_CC[order(rec_omim_CC$x),]
rec_omim_CC$inh <- "rec"
rec_omim_CC$perc <- rec_omim_CC$x/rec_omim_CC$X
omim_CC<-rbind(dom_omim_CC,rec_omim_CC)

vague <- c("intracellular","cell","organelle","cellular_component")
omim_CC <- omim_CC[-which(omim_CC$Description%in%vague),]
omim_CC <- omim_CC[which(omim_CC$corr.p.value<5e-9),]

omim_CC$no <- c(1:length(omim_CC$GO.ID))
omim_CC$no <- factor(omim_CC$no,levels=omim_CC$no)

d<-ggplot(data = omim_CC, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(omim_CC$Description, width = 1, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))


dom_omim_MF <- domomim[which(domomim$category=="MF"),]
dom_omim_MF <- dom_omim_MF[order(dom_omim_MF$x,decreasing=TRUE),]
dom_omim_MF$inh <- "dom"
dom_omim_MF$perc <- dom_omim_MF$x/dom_omim_MF$X

rec_omim_MF <- recomim[which(recomim$category=="MF"),]
rec_omim_MF <- rec_omim_MF[order(rec_omim_MF$x),]
rec_omim_MF$inh <- "rec"
rec_omim_MF$perc <- rec_omim_MF$x/rec_omim_MF$X
omim_MF<-rbind(dom_omim_MF,rec_omim_MF)

vague <- c("molecular_function","binding","protein binding","kinase activity")
omim_MF <- omim_MF[-which(omim_MF$Description%in%vague),]
omim_MF <- omim_MF[which(omim_MF$corr.p.value<5e-9),]

omim_MF$no <- c(1:length(omim_MF$GO.ID))
omim_MF$no <- factor(omim_MF$no,levels=omim_MF$no)

e<-ggplot(data = omim_MF[which(omim_MF$corr.p.value<5e-9),], aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(omim_MF$Description, width = 20, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))



dom_omim_BP <- domomim[which(domomim$category=="BP"),]
dom_omim_BP <- dom_omim_BP[order(dom_omim_BP$x,decreasing=TRUE),]
dom_omim_BP$inh <- "dom"
dom_omim_BP$perc <- dom_omim_BP$x/dom_omim_BP$X

rec_omim_BP <- recomim[which(recomim$category=="BP"),]
rec_omim_BP <- rec_omim_BP[order(rec_omim_BP$x),]
rec_omim_BP$inh <- "rec"
rec_omim_BP$perc <- rec_omim_BP$x/rec_omim_BP$X
omim_BP<-rbind(dom_omim_BP,rec_omim_BP)

vague <- c("biological_process","primary metabolic process","metabolic process","regulation of biological process")
omim_BP <- omim_BP[-which(omim_BP$Description%in%vague),]

omim_BP <- omim_BP[which(omim_BP$corr.p.value<5e-9),]


dups<-omim_BP$Description[which(duplicated(omim_BP$Description)=="TRUE")]
omim_BP <- omim_BP[-which(omim_BP$Description%in%dups),]

omim_BP$no <- c(1:length(omim_BP$GO.ID))
omim_BP$no <- factor(omim_BP$no,levels=omim_BP$no)

f<-ggplot(data = omim_BP, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(dom="#008b45",rec="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(omim_BP$Description, width = 25, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))




ggarrange(a,b,c,d,e,f,ncol=3,nrow=2,widths=c(1,1,1.3))
ggsave("output/Figures/newheatmap.pdf",height=25 , width=16, units='cm')

#CANDIDATES####
write.table(universe_df$gene[which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|
                                                              universe_df$cell_essential=="Y"))],
            "output/genelists/can.txt",row.names = FALSE,col.names=FALSE, quote = FALSE,sep="\t")

candi<-read.table("Gene_lists/BINGO/candidates_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
candi$Genes.in.test.set <- lapply(candi$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

candi$GO.ID <- lapply(candi$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))


candi$category <- vlookup(unlist(candi$GO.ID),catdf,lookup_column = "id",result_column = "cat")
candi$category[which(is.na(candi$category))] <- c("CC","MF","MF")

candi_CC <- candi[which(candi$category=="CC"),]
candi_CC <- candi_CC[order(candi_CC$x),]
candi_CC$perc <- candi_CC$x/candi_CC$X

vague <- c("intracellular","cell","organelle","cellular_component")
candi_CC <- candi_CC[-which(candi_CC$Description%in%vague),]

candi_CC$no <- c(1:length(candi_CC$GO.ID))
candi_CC$no <- factor(candi_CC$no,levels=candi_CC$no)

g<-ggplot(data = candi_CC, aes(x= category,y = no,fill=category,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c("#4594cD"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(candi_CC$Description, width = 1, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))

candi_MF <- candi[which(candi$category=="MF"),]
candi_MF <- candi_MF[order(candi_MF$x),]
candi_MF$perc <- candi_MF$x/candi_MF$X

vague <- c("molecular_function","binding","protein binding","kinase activity")
candi_MF <- candi_MF[-which(candi_MF$Description%in%vague),]

candi_MF$no <- c(1:length(candi_MF$GO.ID))
candi_MF$no <- factor(candi_MF$no,levels=candi_MF$no)

h<-ggplot(data = candi_MF, aes(x= category,y = no,fill=category,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c("#4594cD"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(candi_MF$Description, width = 20, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))

candi_BP <- candi[which(candi$category=="BP"),]
candi_BP <- candi_BP[order(candi_BP$x),]
candi_BP$perc <- candi_BP$x/candi_BP$X

vague <- c("biological_process","primary metabolic process","metabolic process","regulation of biological process")
candi_BP <- candi_BP[-which(candi_BP$Description%in%vague),]

candi_BP$no <- c(1:length(candi_BP$GO.ID))
candi_BP$no <- factor(candi_BP$no,levels=candi_BP$no)

i<-ggplot(data = candi_BP, aes(x= category,y = no,fill=category,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c("#4594cD"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(candi_BP$Description, width = 25, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))


ggarrange(a,b,c,d,e,f,g,h,i,ncol=3,nrow=3,widths=c(1,1.2,1.3))
ggsave("output/Figures/newheatmap_all.pdf",height=22 , width=15, units='cm')
