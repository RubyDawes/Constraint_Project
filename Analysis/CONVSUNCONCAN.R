candi<-read.table("Gene_lists/BINGO/candidates_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
candi$Genes.in.test.set <- lapply(candi$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

candi$GO.ID <- lapply(candi$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))

goslim <- get_ontology("Gene_lists/GOSLIM/go.obo", propagate_relationships = "is_a",extract_tags = "minimal")

category <- unlist(lapply(goslim$ancestors,function(x) ifelse("GO:0003674"%in%x,"MF",
                                                              ifelse("GO:0008150"%in%x,"BP",
                                                                     ifelse("GO:0005575"%in%x,"CC",NA)))))
catdf <- data.frame(id = names(category),cat = unname(category))
rm(category,goslim)

candi$category <- vlookup(unlist(candi$GO.ID),catdf,lookup_column = "id",result_column = "cat")
candi$category[which(is.na(candi$category))] <- c("CC","MF","MF")


candi$no_constrained <- unlist(lapply(1:length(candi$Description),function(x){
  length(which(universe_df$constrained[which(universe_df$gene%in%candi$Genes.in.test.set[[x]])]=="Y"))
}))

candi$no_unconstrained <- unlist(lapply(1:length(candi$Description),function(x){
  length(which(universe_df$constrained[which(universe_df$gene%in%candi$Genes.in.test.set[[x]])]=="N"))
}))


candigenes <- unique(unlist(candi$Genes.in.test.set))

candi$genes_notincat <- lapply(candi$Genes.in.test.set,function(x) candigenes[which(!candigenes%in%x)])
candi$no_notincat_constrained <- unlist(lapply(1:length(candi$Description),function(x){
  length(which(universe_df$constrained[which(universe_df$gene%in%candi$genes_notincat[[x]])]=="Y"))
}))
candi$no_notincat_unconstrained <- unlist(lapply(1:length(candi$Description),function(x){
  length(which(universe_df$constrained[which(universe_df$gene%in%candi$genes_notincat[[x]])]=="N"))
}))

candi$OR <- unlist(lapply(1:length(candi$GO.ID),function(x) {
  fish<-fisher.test(matrix(c(candi$no_constrained[x],candi$no_unconstrained[x],
                             candi$no_notincat_constrained[x],candi$no_notincat_unconstrained[x]
  ),ncol=2),alternative="two.sided")
  return(fish$estimate[[1]])
}))
candi$confintL <- unlist(lapply(1:length(candi$GO.ID),function(x) {
  fish<-fisher.test(matrix(c(candi$no_constrained[x],candi$no_unconstrained[x],
                             candi$no_notincat_constrained[x],candi$no_notincat_unconstrained[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[[1]])
}))
candi$confintH <- unlist(lapply(1:length(candi$GO.ID),function(x) {
  fish<-fisher.test(matrix(c(candi$no_constrained[x],candi$no_unconstrained[x],
                             candi$no_notincat_constrained[x],candi$no_notincat_unconstrained[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[[2]])
}))
candi$pval <- unlist(lapply(1:length(candi$GO.ID),function(x) {
  fish<-fisher.test(matrix(c(candi$no_constrained[x],candi$no_unconstrained[x],
                             candi$no_notincat_constrained[x],candi$no_notincat_unconstrained[x]
  ),ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


min<-min(candi$OR)
max<-max(candi$OR)

candi_CC <- candi[which(candi$category=="CC"),]
candi_CC <- candi_CC[order(candi_CC$OR,decreasing=TRUE),]

vague <- c("intracellular","cell","organelle","cellular_component")
candi_CC <- candi_CC[-which(candi_CC$Description%in%vague),]
candi_CC <- candi_CC[which(candi_CC$pval<5e-2),]

candi_CC$Description<-factor(candi_CC$Description,levels=candi_CC$Description)
candi_CC$gene <- rep("gene",length(candi_CC$GO.ID))

a<-ggplot(data = candi_CC, aes(x= gene,y = Description)) +
  geom_tile(aes(fill = OR,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))


candi_MF <- candi[which(candi$category=="MF"),]
candi_MF <- candi_MF[order(candi_MF$OR,decreasing=TRUE),]

vague <- c("molecular_function","binding","protein binding","kinase activity")
candi_MF <- candi_MF[-which(candi_MF$Description%in%vague),]
candi_MF <- candi_MF[which(candi_MF$pval<5e-2),]

candi_MF$Description<-factor(candi_MF$Description,levels=candi_MF$Description)
candi_MF$gene <- rep("gene",length(candi_MF$GO.ID))

b<-ggplot(data = candi_MF, aes(x= gene,y = Description)) +
  geom_tile(aes(fill = OR,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

candi_BP <- candi[which(candi$category=="BP"),]
candi_BP <- candi_BP[order(candi_BP$OR,decreasing=TRUE),]

vague <- c("biological_process","primary metabolic process","metabolic process","regulation of biological process")
candi_BP <- candi_BP[-which(candi_BP$Description%in%vague),]
candi_BP <- candi_BP[which(candi_BP$pval<5e-2),]

candi_BP$Description<-factor(candi_BP$Description,levels=candi_BP$Description)
candi_BP$gene <- rep("gene",length(candi_BP$GO.ID))

c<-ggplot(data = candi_BP, aes(x= gene,y = Description)) +
  geom_tile(aes(fill = OR,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))





ggarrange(a,b,c,ncol=3,widths=c(1,1,1.5))
ggsave("output/Figures/candidates_strat.pdf",height=12 , width=16, units='cm')



