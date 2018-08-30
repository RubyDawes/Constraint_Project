###############Do lethal genes have different inheritances?#########

lethal_omim<- universe_df[which(universe_df$human_lethal=="Y"&!is.na(universe_df$lethal_inheritance)),c("gene","lethal_inheritance")]
nonlethal_omim<- universe_df[which(universe_df$human_lethal=="N"&universe_df$omim=="Y"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]

lethal_omim$Inheritance_grouped <- ifelse((lethal_omim$lethal_inheritance=="AR"|lethal_omim$lethal_inheritance=="XLr"|lethal_omim$lethal_inheritance=="MT,AR"),"recessive",
                                             ifelse((lethal_omim$lethal_inheritance=="AD"|lethal_omim$lethal_inheritance=="XLd"|lethal_omim$lethal_inheritance=="MT,XLd"),"dominant","recessive/dominant"))
nonlethal_omim$Inheritance_pattern[which(grepl(pattern="MT",nonlethal_omim$Inheritance_pattern))]<-"MT"
nonlethal_omim$Inheritance_pattern[which(grepl("XL",nonlethal_omim$Inheritance_pattern)&!grepl("MT",nonlethal_omim$Inheritance_pattern))]<-"XL"
nonlethal_omim$Inheritance_grouped <- ifelse((nonlethal_omim$Inheritance_pattern=="AR"|nonlethal_omim$Inheritance_pattern=="XLr"|nonlethal_omim$Inheritance_pattern=="MT"),"recessive",
                                               ifelse((nonlethal_omim$Inheritance_pattern=="AD"|nonlethal_omim$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))


lethal_omim_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                         value=c(length(which(lethal_omim$Inheritance_grouped=="dominant")),
                                                 length(which(lethal_omim$Inheritance_grouped=="recessive")),
                                                 length(which(lethal_omim$Inheritance_grouped=="recessive/dominant"))))
lethal_omim_inh_grouped$group <- factor(lethal_omim_inh_grouped$group, levels = lethal_omim_inh_grouped$group)
a <- ggplot(lethal_omim_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("human lethal genes"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))

nonlethal_omim_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                           value=c(length(which(nonlethal_omim$Inheritance_grouped=="dominant")),
                                                   length(which(nonlethal_omim$Inheritance_grouped=="recessive")),
                                                   length(which(nonlethal_omim$Inheritance_grouped=="recessive/dominant"))))
nonlethal_omim_inh_grouped$group <- factor(nonlethal_omim_inh_grouped$group, levels = nonlethal_omim_inh_grouped$group)
b <- ggplot(nonlethal_omim_inh_grouped, aes(x="", y=value, fill=group))
b<- b+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
b<- b+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
b<- b+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
b<- b+ggtitle(paste("Non-lethal OMIM genes"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))


ggsave("Analysis/Sandra_Figures/Figs/inheritance_lethalvsnonlethal.pdf",arrangeGrob(a,b,ncol=2),height=9, width=16, units='cm')
