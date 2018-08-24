#####################1B: inheritance stacked bar charts####################
const <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$constrained=="Y"))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$constrained=="N"))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))),
                      AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))),
                      "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))),
                      AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AD"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))),
                      XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$constrained=="Y"))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$constrained=="N"))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))))

constm <- melt(const)
constm$constraint <- rep(c("Constraint     ", "No constraint     "), 5)

f <- ggplot(dat=constm, aes(x=variable, y=value, fill=constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("#CD5555","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/whole_gene_constraint.pdf",height=7, width=7, units='cm')
rm(const,constm,f)
regconst <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$any_constraint))),
                           length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$any_constraint)))),
                      AR=c(length(which(universe_df$Inheritance_pattern=="AR"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$any_constraint))),
                           length(which(universe_df$Inheritance_pattern=="AR"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$any_constraint)))),
                      "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$any_constraint))),
                               length(which(universe_df$Inheritance_pattern=="AR,AD"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$any_constraint)))),
                      AD=c(length(which(universe_df$Inheritance_pattern=="AD"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$any_constraint))),
                           length(which(universe_df$Inheritance_pattern=="AD"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$any_constraint)))),
                      XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$any_constraint))),
                           length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$any_constraint)))))


regconstm <- melt(regconst)
regconstm$regconstraint <- rep(c("Regional or whole-gene constraint     ", "no regional or whole-gene constraint     "), 5)

e <- ggplot(dat=regconstm, aes(x=variable, y=value, fill=regconstraint))
e<- e+geom_bar(width = 0.8, stat = "identity",color="slategray")+scale_y_continuous(expand = c(0, 0)) +bar_theme()
e<- e+labs(y = "Proportion")+scale_fill_manual(values=c("steelblue3","lightsalmon2"))+theme(legend.position="bottom")


ggsave("Analysis/Sandra_Figures/Figs/reg_constraint.pdf",height=7, width=7, units='cm')
rm(regconst,regconstm,e)


#################pie of inheritance of constrained vs non-constrained genes###############

constrained<- universe_df[which(universe_df$constrained=="Y"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]
unconstrained<- universe_df[which(universe_df$constrained=="N"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]

constrained$Inheritance_pattern[which(grepl(pattern="MT",constrained$Inheritance_pattern))]<-"MT"
constrained$Inheritance_pattern[which(grepl("XL",constrained$Inheritance_pattern)&!grepl("MT",constrained$Inheritance_pattern))]<-"XL"
constrained$Inheritance_grouped <- ifelse((constrained$Inheritance_pattern=="AR"|constrained$Inheritance_pattern=="XLr"|constrained$Inheritance_pattern=="MT"),"recessive",
                                          ifelse((constrained$Inheritance_pattern=="AD"|constrained$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))
unconstrained$Inheritance_pattern[which(grepl(pattern="MT",unconstrained$Inheritance_pattern))]<-"MT"
unconstrained$Inheritance_pattern[which(grepl("XL",unconstrained$Inheritance_pattern)&!grepl("MT",unconstrained$Inheritance_pattern))]<-"XL"
unconstrained$Inheritance_grouped <- ifelse((unconstrained$Inheritance_pattern=="AR"|unconstrained$Inheritance_pattern=="XLr"|unconstrained$Inheritance_pattern=="MT"),"recessive",
                                          ifelse((unconstrained$Inheritance_pattern=="AD"|unconstrained$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))


constrained_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                      value=c(length(which(constrained$Inheritance_grouped=="dominant")),
                                              length(which(constrained$Inheritance_grouped=="recessive")),
                                              length(which(constrained$Inheritance_grouped=="recessive/dominant"))))
constrained_inh_grouped$group <- factor(constrained_inh_grouped$group, levels = constrained_inh_grouped$group)
a <- ggplot(constrained_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Genes with whole-gene constraint"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))
ggsave("Analysis/Sandra_Figures/Figs/constrained_inheritance_pie.pdf",height=7, width=7, units='cm')

unconstrained_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                      value=c(length(which(unconstrained$Inheritance_grouped=="dominant")),
                                              length(which(unconstrained$Inheritance_grouped=="recessive")),
                                              length(which(unconstrained$Inheritance_grouped=="recessive/dominant"))))
unconstrained_inh_grouped$group <- factor(unconstrained_inh_grouped$group, levels = unconstrained_inh_grouped$group)
a <- ggplot(unconstrained_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Genes with no whole-gene constraint"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))
ggsave("Analysis/Sandra_Figures/Figs/unconstrained_inheritance_pie.pdf",height=7, width=7, units='cm')


#####################doing the same but including regional constraint##################
constrained<- universe_df[which(universe_df$any_constraint=="Y"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]
unconstrained<- universe_df[which(universe_df$any_constraint=="N"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]

constrained$Inheritance_pattern[which(grepl(pattern="MT",constrained$Inheritance_pattern))]<-"MT"
constrained$Inheritance_pattern[which(grepl("XL",constrained$Inheritance_pattern)&!grepl("MT",constrained$Inheritance_pattern))]<-"XL"
constrained$Inheritance_grouped <- ifelse((constrained$Inheritance_pattern=="AR"|constrained$Inheritance_pattern=="XLr"|constrained$Inheritance_pattern=="MT"),"recessive",
                                          ifelse((constrained$Inheritance_pattern=="AD"|constrained$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))
unconstrained$Inheritance_pattern[which(grepl(pattern="MT",unconstrained$Inheritance_pattern))]<-"MT"
unconstrained$Inheritance_pattern[which(grepl("XL",unconstrained$Inheritance_pattern)&!grepl("MT",unconstrained$Inheritance_pattern))]<-"XL"
unconstrained$Inheritance_grouped <- ifelse((unconstrained$Inheritance_pattern=="AR"|unconstrained$Inheritance_pattern=="XLr"|unconstrained$Inheritance_pattern=="MT"),"recessive",
                                            ifelse((unconstrained$Inheritance_pattern=="AD"|unconstrained$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))


constrained_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                      value=c(length(which(constrained$Inheritance_grouped=="dominant")),
                                              length(which(constrained$Inheritance_grouped=="recessive")),
                                              length(which(constrained$Inheritance_grouped=="recessive/dominant"))))
constrained_inh_grouped$group <- factor(constrained_inh_grouped$group, levels = constrained_inh_grouped$group)
a <- ggplot(constrained_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Genes with regional or whole-gene constraint"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))
ggsave("Analysis/Sandra_Figures/Figs/reg_and_wholeconstrained_inheritance_pie.pdf",height=7, width=7, units='cm')

unconstrained_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                        value=c(length(which(unconstrained$Inheritance_grouped=="dominant")),
                                                length(which(unconstrained$Inheritance_grouped=="recessive")),
                                                length(which(unconstrained$Inheritance_grouped=="recessive/dominant"))))
unconstrained_inh_grouped$group <- factor(unconstrained_inh_grouped$group, levels = unconstrained_inh_grouped$group)
a <- ggplot(unconstrained_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Genes with no regional or whole-gene constraint"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))
ggsave("Analysis/Sandra_Figures/Figs/reg_and_wholeunconstrained_inheritance_pie.pdf",height=7, width=7, units='cm')


###############pie of inheritance of constrained vs non-constrained HUMAN LETHAL genes#########

constrained<- universe_df[which(universe_df$constrained=="Y"&universe_df$human_lethal=="Y"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]
unconstrained<- universe_df[which(universe_df$constrained=="N"&universe_df$human_lethal=="Y"&!is.na(universe_df$Inheritance_pattern)),c("gene","Inheritance_pattern")]

constrained$Inheritance_pattern[which(grepl(pattern="MT",constrained$Inheritance_pattern))]<-"MT"
constrained$Inheritance_pattern[which(grepl("XL",constrained$Inheritance_pattern)&!grepl("MT",constrained$Inheritance_pattern))]<-"XL"
constrained$Inheritance_grouped <- ifelse((constrained$Inheritance_pattern=="AR"|constrained$Inheritance_pattern=="XLr"|constrained$Inheritance_pattern=="MT"),"recessive",
                                          ifelse((constrained$Inheritance_pattern=="AD"|constrained$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))
unconstrained$Inheritance_pattern[which(grepl(pattern="MT",unconstrained$Inheritance_pattern))]<-"MT"
unconstrained$Inheritance_pattern[which(grepl("XL",unconstrained$Inheritance_pattern)&!grepl("MT",unconstrained$Inheritance_pattern))]<-"XL"
unconstrained$Inheritance_grouped <- ifelse((unconstrained$Inheritance_pattern=="AR"|unconstrained$Inheritance_pattern=="XLr"|unconstrained$Inheritance_pattern=="MT"),"recessive",
                                            ifelse((unconstrained$Inheritance_pattern=="AD"|unconstrained$Inheritance_pattern=="XLd"),"dominant","recessive/dominant"))


constrained_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                      value=c(length(which(constrained$Inheritance_grouped=="dominant")),
                                              length(which(constrained$Inheritance_grouped=="recessive")),
                                              length(which(constrained$Inheritance_grouped=="recessive/dominant"))))
constrained_inh_grouped$group <- factor(constrained_inh_grouped$group, levels = constrained_inh_grouped$group)
a <- ggplot(constrained_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Prenatal/infantile lethal genes with depleted genetic variation"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))
ggsave("Analysis/Sandra_Figures/Figs/lethal_constrained_inheritance_pie.pdf",height=7, width=7, units='cm')

unconstrained_inh_grouped <- data.frame(group=c("Dominant     ","Recessive     ","Dominant/Recessive"),
                                        value=c(length(which(unconstrained$Inheritance_grouped=="dominant")),
                                                length(which(unconstrained$Inheritance_grouped=="recessive")),
                                                length(which(unconstrained$Inheritance_grouped=="recessive/dominant"))))
unconstrained_inh_grouped$group <- factor(unconstrained_inh_grouped$group, levels = unconstrained_inh_grouped$group)
a <- ggplot(unconstrained_inh_grouped, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_brewer(palette="Accent")+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Prenatal/infantile lethal genes with normal genetic variation"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))
ggsave("Analysis/Sandra_Figures/Figs/lethal_unconstrained_inheritance_pie.pdf",height=7, width=7, units='cm')







