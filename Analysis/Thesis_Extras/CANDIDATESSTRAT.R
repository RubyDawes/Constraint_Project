#candidates####

write.table(universe_df$gene[which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|
                                                              universe_df$cell_essential=="Y")&universe_df$constrained=="Y")],
            "output/genelists/con_can.txt",row.names = FALSE,col.names=FALSE, quote = FALSE,sep="\t")
write.table(universe_df$gene[which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|
                                                              universe_df$cell_essential=="Y")&universe_df$constrained=="N")],
            "output/genelists/uncon_can.txt",row.names = FALSE,col.names=FALSE, quote = FALSE,sep="\t")


uncon_can<-read.table("Gene_lists/BINGO/unconstrained_candidates_0.05.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
uncon_can$Genes.in.test.set <- lapply(uncon_can$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

con_can<-read.table("Gene_lists/BINGO/constrained_candidates_0.05.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
con_can$Genes.in.test.set <- lapply(con_can$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

uncon_can$GO.ID <- lapply(uncon_can$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))
con_can$GO.ID <- lapply(con_can$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))


uncon_can$category <- vlookup(unlist(uncon_can$GO.ID),catdf,lookup_column = "id",result_column = "cat")
uncon_can$category[which(is.na(uncon_can$category))] <- c("CC","MF","MF")
con_can$category <- vlookup(unlist(con_can$GO.ID),catdf,lookup_column = "id",result_column = "cat")
con_can$category[which(is.na(con_can$category))] <- c("CC","MF","MF","CC")

con_can_CC <- con_can[which(con_can$category=="CC"),]
con_can_CC <- con_can_CC[order(con_can_CC$x,decreasing=TRUE),]
con_can_CC$inh <- "con"
con_can_CC$perc <- con_can_CC$x/con_can_CC$X

uncon_can_CC <- uncon_can[which(uncon_can$category=="CC"),]
uncon_can_CC <- uncon_can_CC[order(uncon_can_CC$x),]
uncon_can_CC$inh <- "uncon"
uncon_can_CC$perc <- uncon_can_CC$x/uncon_can_CC$X
candidates_CC<-rbind(con_can_CC,uncon_can_CC)

vague <- c("intracellular","cell","organelle","cellular_component")
candidates_CC <- candidates_CC[-which(candidates_CC$Description%in%vague),]
candidates_CC <- candidates_CC[which(candidates_CC$corr.p.value<5e-9),]

candidates_CC$no <- c(1:length(candidates_CC$GO.ID))
candidates_CC$no <- factor(candidates_CC$no,levels=candidates_CC$no)

d<-ggplot(data = candidates_CC, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(con="#008b45",uncon="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(candidates_CC$Description, width = 1, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))


con_can_MF <- con_can[which(con_can$category=="MF"),]
con_can_MF <- con_can_MF[order(con_can_MF$x,decreasing=TRUE),]
con_can_MF$inh <- "con"
con_can_MF$perc <- con_can_MF$x/con_can_MF$X

uncon_can_MF <- uncon_can[which(uncon_can$category=="MF"),]
uncon_can_MF <- uncon_can_MF[order(uncon_can_MF$x),]
uncon_can_MF$inh <- "uncon"
uncon_can_MF$perc <- uncon_can_MF$x/uncon_can_MF$X
candidates_MF<-rbind(con_can_MF,uncon_can_MF)

vague <- c("molecular_function","binding","protein binding","kinase activity")
candidates_MF <- candidates_MF[-which(candidates_MF$Description%in%vague),]
candidates_MF <- candidates_MF[which(candidates_MF$corr.p.value<5e-9),]

candidates_MF$no <- c(1:length(candidates_MF$GO.ID))
candidates_MF$no <- factor(candidates_MF$no,levels=candidates_MF$no)

d<-ggplot(data = candidates_MF, aes(x= category,y = no,fill=inh,alpha=perc)) +
  geom_tile(aes(width=1)) +
  scale_fill_manual(values=c(con="#008b45",uncon="#a020f0"))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position="none",axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(candidates_MF$Description, width = 1, simplify = FALSE), paste, collapse="\n"))+
  scale_x_discrete(expand = c(0, 0))+scale_alpha(limits = c(0,0.7))




