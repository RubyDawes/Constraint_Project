############pie of human lethal genes##############
human_lethal_pie <- data.frame(group=c("human Lethal","OMIM Non-Lethal", "mouse lethal non-OMIM","non-OMIM,non-lethal"),
                               value=c(length(which(universe_df$human_lethal=="Y")),
                                       length(which(universe_df$omim=="Y"&universe_df$human_lethal=="N")),
                                       length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                                       length(which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="N"|is.na(universe_df$lethal_mouse))))
                                       
                                       ))
human_lethal_pie$group <- factor(human_lethal_pie$group, levels = human_lethal_pie$group)
c <- ggplot(human_lethal_pie, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("#961b28","#d36e7d","black","#52b7d2"))
c<- c+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
ggsave("Analysis/Sandra_Figures/Figs/humanlethal_pie_perc.pdf",height=18, width=18, units='cm')

##############pie for hospital talk################
human_lethal_pie <- data.frame(group=c("Disease genes","Lethal","non-Disease"),
                               value=c(length(which(universe_df$omim=="Y"&universe_df$human_lethal=="N")),
                                       length(which(universe_df$human_lethal=="Y")),
                                       length(which(is.na(universe_df$omim)))))
human_lethal_pie$group <- factor(human_lethal_pie$group, levels = human_lethal_pie$group)
c <- ggplot(human_lethal_pie, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("#d36e7d","#961b28","steelblue3"))
c<- c+geom_text(aes(y = value,label = paste(percent(value/sum(value)),"\n",value,sep="")),size=6,position = position_stack(vjust = 0.5))
ggsave("Analysis/Sandra_Figures/Figs/hospweek_pie_one.pdf",height=18, width=18, units='cm')

human_lethal_pie <- data.frame(group=c("Disease genes","Lethal","Mouse Lethal, No known phenotype in humans","non-Disease"),
                               value=c(length(which(universe_df$omim=="Y"&universe_df$human_lethal_B=="N")),
                                       length(which(universe_df$human_lethal_B=="Y")),
                                       length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                                       length(which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="N"|is.na(universe_df$lethal_mouse))))))
human_lethal_pie$group <- factor(human_lethal_pie$group, levels = human_lethal_pie$group)
c <- ggplot(human_lethal_pie, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("#d36e7d","#961b28","black","steelblue3"))
c<- c+geom_text(aes(y = value,label = paste(percent(value/sum(value)),"\n",value,sep="")),size=6,position = position_stack(vjust = 0.5))
ggsave("Analysis/Sandra_Figures/Figs/hospweek_pie_two.pdf",height=18, width=18, units='cm')

############pie of human lethal genes##############
human_lethal_pie <- data.frame(group=c("omim", "non-omim"),
                               value=c(
                                       length(which(universe_df$omim=="Y")),
                                       length(which(is.na(universe_df$omim)))))
human_lethal_pie$group <- factor(human_lethal_pie$group, levels = human_lethal_pie$group)
c <- ggplot(human_lethal_pie, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("#961b28","#d36e7d","black","#52b7d2"))
c<- c+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
ggsave("Analysis/Sandra_Figures/Figs/midyearpie.pdf",height=18, width=18, units='cm')

##########Figure 3: Dominant and X-linked lethal genes are significantly more intolerant to genetic variation,than recessive genes##########
lethal_omim<- universe_df[which(universe_df$human_lethal=="Y"&!is.na(universe_df$lethal_inheritance)),c("gene","lethal_inheritance","constrained")]
lethal_omim<-lethal_omim[-which(is.na(lethal_omim$constrained)),]
lethal_omim$Inheritance_grouped <- ifelse((lethal_omim$lethal_inheritance=="AR"|lethal_omim$lethal_inheritance=="XLr"|lethal_omim$lethal_inheritance=="MT,AR"),"recessive",
                                          ifelse((lethal_omim$lethal_inheritance=="AD"|lethal_omim$lethal_inheritance=="XLd"|lethal_omim$lethal_inheritance=="MT,XLd"),"dominant","recessive/dominant"))

lethal_constraint <- data.frame(Recessive=c(length(which(lethal_omim$Inheritance_grouped=="recessive"&lethal_omim$constrained=="Y"))/length(which(lethal_omim$Inheritance_grouped=="recessive")),
                                length(which(lethal_omim$Inheritance_grouped=="recessive"&lethal_omim$constrained=="N"))/length(which(lethal_omim$Inheritance_grouped=="recessive"))),
                                Mixed=c(length(which(lethal_omim$Inheritance_grouped=="recessive/dominant"&lethal_omim$constrained=="Y"))/length(which(lethal_omim$Inheritance_grouped=="recessive/dominant")),
                               length(which(lethal_omim$Inheritance_grouped=="recessive/dominant"&lethal_omim$constrained=="N"))/length(which(lethal_omim$Inheritance_grouped=="recessive/dominant"))),
                               Dominant=c(length(which(lethal_omim$Inheritance_grouped=="dominant"&lethal_omim$constrained=="Y"))/length(which(lethal_omim$Inheritance_grouped=="dominant")),
                                          length(which(lethal_omim$Inheritance_grouped=="dominant"&lethal_omim$constrained=="N"))/length(which(lethal_omim$Inheritance_grouped=="dominant"))))

lethal_constraintm <- melt(lethal_constraint)
lethal_constraintm$constraint <- rep(c("Constraint     ", "No constraint     "), 3)

f <- ggplot(dat=lethal_constraintm, aes(x=variable, y=value, fill=constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("#CD5555","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/grouped_constraint.pdf",height=7, width=7, units='cm')
