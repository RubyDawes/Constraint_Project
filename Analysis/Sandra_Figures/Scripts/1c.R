
############yeast##########
yeast_viable <-read.table("Gene_lists/EG/yeast/viable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable <-read.table("Gene_lists/EG/yeast/inviable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable<-yeast_inviable[which(yeast_inviable$V6=="null"),]
yeast_viable<-yeast_viable[which(yeast_viable$V6=="null"),]

lethal_prop <- data.frame(Yeast=c(length(unique(yeast_inviable$V1))/(length(unique(yeast_inviable$V1))+length(unique(yeast_viable$V1))),
                                  length(unique(yeast_viable$V1))/(length(unique(yeast_inviable$V1))+length(unique(yeast_viable$V1)))),
                                Cells=c(length(which(universe_df$cell_essential_hits>=2))/length(which(!is.na(universe_df$cell_ko))),
                                        length(which(universe_df$cell_essential_hits<2))/length(which(!is.na(universe_df$cell_ko)))),
                                Mouse=c(length(which(universe_df$lethal_mouse=="Y"))/length(which(universe_df$mouse_ko=="Y")),
                                        length(which(universe_df$lethal_mouse=="N"))/length(which(universe_df$mouse_ko=="Y"))),
                                Human=c(length(which(universe_df$human_lethal_B=="Y"))/length(universe_df$gene),
                                        length(which(universe_df$human_lethal_B=="N"))/length(universe_df$gene)))

lethal_propm <- melt(lethal_prop)
lethal_propm$mis_constraint <- rep(c("Lethal     ","Non-lethal     "), 4)

f <- ggplot(dat=lethal_propm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="right")
ggsave("Analysis/Sandra_Figures/Figs/lethal_prop_1c.pdf",height=9, width=20.5, units='cm')
