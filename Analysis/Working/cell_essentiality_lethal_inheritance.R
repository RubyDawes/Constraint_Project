###############Are recessive disease genes more likely to be cell essential?#########

cell_essential<- universe_df[which(universe_df$human_lethal=="Y"&universe_df$cell_essential=="Y"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]
uncell_essential<- universe_df[which(universe_df$human_lethal=="Y"&universe_df$cell_essential=="N"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]

cell_essential$Inheritance_pattern[which(grepl(pattern="MT",cell_essential$Inheritance_pattern))]<-"MT"
cell_essential$Inheritance_pattern[which(grepl("XL",cell_essential$Inheritance_pattern)&!grepl("MT",cell_essential$Inheritance_pattern))]<-"XL"
cell_essential$Inheritance_grouped <- ifelse((cell_essential$Inheritance_pattern=="AR"|cell_essential$Inheritance_pattern=="XLr"|cell_essential$Inheritance_pattern=="MT"),"recessive",
                                          ifelse((cell_essential$Inheritance_pattern=="AD"|cell_essential$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))
uncell_essential$Inheritance_pattern[which(grepl(pattern="MT",uncell_essential$Inheritance_pattern))]<-"MT"
uncell_essential$Inheritance_pattern[which(grepl("XL",uncell_essential$Inheritance_pattern)&!grepl("MT",uncell_essential$Inheritance_pattern))]<-"XL"
uncell_essential$Inheritance_grouped <- ifelse((uncell_essential$Inheritance_pattern=="AR"|uncell_essential$Inheritance_pattern=="XLr"|uncell_essential$Inheritance_pattern=="MT"),"recessive",
                                            ifelse((uncell_essential$Inheritance_pattern=="AD"|uncell_essential$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))


cell_essential_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                      value=c(length(which(cell_essential$Inheritance_grouped=="dominant")),
                                              length(which(cell_essential$Inheritance_grouped=="recessive")),
                                              length(which(cell_essential$Inheritance_grouped=="recessive/dominant"))))
cell_essential_inh_grouped$group <- factor(cell_essential_inh_grouped$group, levels = cell_essential_inh_grouped$group)
a <- ggplot(cell_essential_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Cell essential human lethal genes"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))

uncell_essential_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                         value=c(length(which(uncell_essential$Inheritance_grouped=="dominant")),
                                                 length(which(uncell_essential$Inheritance_grouped=="recessive")),
                                                 length(which(uncell_essential$Inheritance_grouped=="recessive/dominant"))))
uncell_essential_inh_grouped$group <- factor(uncell_essential_inh_grouped$group, levels = uncell_essential_inh_grouped$group)
a <- ggplot(uncell_essential_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Cell unessential human lethal genes"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))



ggsave("Analysis/Sandra_Figures/Figs/lethal_cell_essential_inheritance_pie.pdf",height=7, width=7, units='cm')
