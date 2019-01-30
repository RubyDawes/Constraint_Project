source("plot_functions.R")

#2A: Scatter plot showing levels of genetic tolerance to LoF (pLI) or missense constraint####
nonomim <- ggplot(universe_df[which(is.na(universe_df$omim)&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,alpha=0.3)+scatter_theme()+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("Non-OMIM genes \n n=",length(which(is.na(universe_df$omim)&!is.na(universe_df$exac))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01))+
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14))

omim <- ggplot(universe_df[which(universe_df$omim=="Y"&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,alpha=0.3)+scatter_theme()+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("OMIM genes \n n=",length(which(universe_df$omim=="Y"&!is.na(universe_df$exac))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14)) 

ggarrange(omim,nonomim)
ggsave("output/Figures/2A.pdf",height=7.6, width=17.5, units='cm')

rm(nonomim,omim)

#2B: Pie Charts contrasting relative levels of genetic constraint for OMIM versus non-OMIM genes####

omim_mis <- data.frame(group=c("Missense constraint     ","No missense constraint     "),value=c(length(which(universe_df$omim=="Y"&universe_df$mis_z>=3.09)),
                                                                                                 length(which(universe_df$omim=="Y"&universe_df$mis_z<3.09))))
nonomim_mis <- data.frame(group=c("Missense constraint     ","No missense constraint     "),value=c(length(which(is.na(universe_df$omim)&universe_df$mis_z>=3.09)),
                                                                                                    length(which(is.na(universe_df$omim)&universe_df$mis_z<3.09))))
omim_mis$group <- factor(omim_mis$group, levels = omim_mis$group)
nonomim_mis$group <- factor(nonomim_mis$group, levels = nonomim_mis$group)
a <- ggplot(omim_mis, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = percent(value/sum(value))),size = 4,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("OMIM genes"))+theme(plot.title = element_text(size = 10, face = "bold"),legend.position="none",legend.text=element_text(size = 10))

b <- ggplot(nonomim_mis, aes(x="", y=value, fill=group))
b<- b+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
b<- b+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="none")
b<- b+geom_text(aes(y = value,label = percent(value/sum(value))),size = 4,position = position_stack(vjust = 0.5))
b<- b+ggtitle(paste("Non-OMIM genes"))+theme(plot.title = element_text(size = 10, face = "bold"),legend.position="none",legend.text=element_text(size = 10))


omim_lof <- data.frame(group=c("LoF constraint     ","No LoF constraint     "),value=c(length(which(universe_df$omim=="Y"&universe_df$pLI>=0.9)),
                                                                                       length(which(universe_df$omim=="Y"&universe_df$pLI<0.9))))
nonomim_lof <- data.frame(group=c("LoF constraint     ","No LoF constraint     "),value=c(length(which(is.na(universe_df$omim)&universe_df$pLI>=0.9)),
                                                                                          length(which(is.na(universe_df$omim)&universe_df$pLI<0.9))))
omim_lof$group <- factor(omim_lof$group, levels = omim_lof$group)
nonomim_lof$group <- factor(nonomim_lof$group, levels = nonomim_lof$group)
c <- ggplot(omim_lof, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="none")
c<- c+geom_text(aes(y = value,label = percent(value/sum(value))),size = 4,position = position_stack(vjust = 0.5))
c<- c+ggtitle(paste("OMIM genes"))+theme(plot.title = element_text(size = 10, face = "bold"),legend.position="none",legend.text=element_text(size = 10))

d <- ggplot(nonomim_lof, aes(x="", y=value, fill=group))
d<- d+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
d<- d+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="none")
d<- d+geom_text(aes(y = value,label = percent(value/sum(value))),size = 4,position = position_stack(vjust = 0.5))
d<- d+ggtitle(paste("Non-OMIM genes"))+theme(plot.title = element_text(size = 10, face = "bold"),legend.position="none",legend.text=element_text(size = 10))

ggarrange(
  ggarrange(a,b,common.legend=TRUE,legend="bottom"),
  ggarrange(c,d,common.legend=TRUE,legend="bottom")
)
ggsave("output/Figures/2B.pdf",height=5, width=21, units='cm')
dev.off()
rm(a,b,c,d,omim_mis,nonomim_mis,omim_lof,nonomim_lof)


#2C: Histograms showing correlation of inheritance pattern with levels of genetic constraint####
levels<- c(paste0("MT \n n = ",length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))),
           paste0("AR \n n = ",length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))),
           paste0("AR/AD \n n = ",length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))),
           paste0("AD \n n = ",length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))),
           paste0("XL \n n = ",length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))))

inh_mis <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&(universe_df$mis_z>=3.09)))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y"))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))),
                      AR=c(length(which(universe_df$Inheritance_pattern=="AR"&(universe_df$mis_z>=3.09)))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AR"&(universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AR"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))),
                      "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&(universe_df$mis_z>=3.09)))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AR,AD"&(universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AR,AD"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))),
                      AD=c(length(which(universe_df$Inheritance_pattern=="AD"&(universe_df$mis_z>=3.09)))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AD"&(universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AD"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))),
                      XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&(universe_df$mis_z>=3.09)))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&(universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y")))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))))
inh_mism <- melt(inh_mis)
inh_mism$mis_constraint <- rep(c("Whole-gene missense constraint     ","Regional missense constraint     ", "No missense constraint     "), 5)
inh_mism$mis_constraint <-factor(inh_mism$mis_constraint,levels=c("Whole-gene missense constraint     ","Regional missense constraint     ", "No missense constraint     "))

inh_mism$variable <- rep(levels,each=3)
inh_mism$variable <- factor(inh_mism$variable, levels = levels)


labelsmis <- rep("",15)
labelsmis[c(2,4,5,7,8,10,11,13,14)]<-c(paste0(round(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(universe_df$Inheritance_pattern=="AR"&(universe_df$mis_z>=3.09)))/
                                                     length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(universe_df$Inheritance_pattern=="AR"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                     length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(universe_df$Inheritance_pattern=="AR,AD"&(universe_df$mis_z>=3.09)))/
                                        length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(universe_df$Inheritance_pattern=="AR,AD"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                     length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(universe_df$Inheritance_pattern=="AD"&(universe_df$mis_z>=3.09)))/
                                        length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(universe_df$Inheritance_pattern=="AD"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                     length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&(universe_df$mis_z>=3.09)))/
                                        length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))*100,0),"%"),
                                      paste0(round(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                        length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))*100,0),"%"))

f <- ggplot(dat=inh_mism, aes(x=variable, y=value, fill=mis_constraint))+
  geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse=TRUE))+scale_y_continuous(expand = c(0, 0),labels = percent) +bar_theme()+
  labs(y = "Percent constrained")+scale_fill_manual(values=c("lightsalmon2","#f9b8a0","steelblue3"))+theme(legend.position="none")+
  geom_text(aes(label = labelsmis),size=3,position=position_fill(reverse=TRUE,vjust=1))


inh_lof <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$pLI>=0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$pLI<0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))),
                      AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))),
                      "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))),
                      AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))),
                      XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$pLI>=0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$pLI<0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))))
inh_lofm <- melt(inh_lof)
inh_lofm$lof_constraint <- rep(c("LoF constraint     ", "no LoF constraint     "), 5)

inh_lofm$variable <- rep(levels,each=2)
inh_lofm$variable <- factor(inh_lofm$variable, levels = levels)
labels <- rep("",10)
labels[c(1,3,5,7,9)]<-c(paste0(round(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$pLI>=0.9))/
                        length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9))/
                          length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI>=0.9))/
                          length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))/
                          length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$pLI>=0.9))/
                          length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))*100,0),"%"))

e <- ggplot(dat=inh_lofm, aes(x=variable, y=value, fill=lof_constraint))
e<- e+geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0),labels = percent) +bar_theme()
e<- e+labs(y = "Percent constrained")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="none")+geom_text(aes(label = labels),size=3,position = "identity",vjust=-1)


ggarrange(f,e)
ggsave("output/Figures/2C.pdf",height=5, width=13, units='cm')


#2D: human lethal genes constraint ####

levels<- c(paste0("MT \n n = ",length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
           paste0("AR \n n = ",length(which((universe_df$lethal_inheritance=="AR")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
           paste0("AR/AD \n n = ",length(which((universe_df$lethal_inheritance=="AR,AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
           paste0("AD \n n = ",length(which((universe_df$lethal_inheritance=="AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
           paste0("XL \n n = ",length(which((grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))))


inh_mis_human_lethal_B <- data.frame(MT=c(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z>=3.09))/length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y"))/length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))/length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AR=c(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR")&universe_df$mis_z>=3.09))/length(which((universe_df$lethal_inheritance=="AR")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y"))/length(which((universe_df$lethal_inheritance=="AR")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))/length(which((universe_df$lethal_inheritance=="AR")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     "AR,AD"=c(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR,AD")&universe_df$mis_z>=3.09))/length(which((universe_df$lethal_inheritance=="AR,AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                              length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR,AD")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y"))/length(which((universe_df$lethal_inheritance=="AR,AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                              length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR,AD")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))/length(which((universe_df$lethal_inheritance=="AR,AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AD=c(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD")&universe_df$mis_z>=3.09))/length(which((universe_df$lethal_inheritance=="AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y"))/length(which((universe_df$lethal_inheritance=="AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD")&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))/length(which((universe_df$lethal_inheritance=="AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     XL=c(length(which(universe_df$human_lethal_B=="Y"&(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$mis_z>=3.09))/length(which((grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="Y"))/length(which((grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$human_lethal_B=="Y"&(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))/length(which((grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))))
                                       
                                       
inh_mis_human_lethal_Bm <- melt(inh_mis_human_lethal_B)
inh_mis_human_lethal_Bm$mis_constraint <- rep(c("Missense constraint     ","Regional missense constraint", "No missense constraint     "), 5)

inh_mis_human_lethal_Bm$variable<- rep(levels,each=3)
inh_mis_human_lethal_Bm$variable <- factor(inh_mis_human_lethal_Bm$variable, levels = levels)

labelsf <- rep("",15)
labelsf[c(2,4,5,7,8,10,11,13,14)]<-c(paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                    length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR")&universe_df$mis_z>=3.09))/
                                                      length(which((universe_df$lethal_inheritance=="AR")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR")&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                      length(which((universe_df$lethal_inheritance=="AR")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR,AD")&universe_df$mis_z>=3.09))/
                                                      length(which((universe_df$lethal_inheritance=="AR,AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR,AD")&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                      length(which((universe_df$lethal_inheritance=="AR,AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD")&universe_df$mis_z>=3.09))/
                                                      length(which((universe_df$lethal_inheritance=="AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD")&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                      length(which((universe_df$lethal_inheritance=="AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$mis_z>=3.09))/
                                                      length(which((grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/
                                                      length(which((grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"))

inh_mis_human_lethal_Bm$mis_constraint <-factor(inh_mis_human_lethal_Bm$mis_constraint,levels=c("Missense constraint     ","Regional missense constraint", "No missense constraint     "))

m <- ggplot(dat=inh_mis_human_lethal_Bm, aes(x=variable, y=value, fill=mis_constraint))+
  geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse = TRUE))+
  scale_y_continuous(expand = c(0, 0),labels=percent) +bar_theme()+
  labs(y = "Percent constrained")+scale_fill_manual(values=c("lightsalmon2","#f9b8a0","steelblue3"))+
  theme(legend.position="none")+geom_text(aes(label = labelsf),size=3,position = position_fill(reverse=TRUE,vjust=1))



inh_lof_human_lethal_B <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AR=c(length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     "ARAD"=c(length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                              length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AD=c(length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     XL=c(length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))))


inh_lof_human_lethal_Bm <- melt(inh_lof_human_lethal_B)
inh_lof_human_lethal_Bm$mis_constraint <- rep(c("LoF constraint     ", "No LoF constraint     "), 5)

inh_lof_human_lethal_Bm$variable<- rep(levels,each=2)
inh_lof_human_lethal_Bm$variable <- factor(inh_lof_human_lethal_Bm$variable, levels = levels)

labels <- rep("",10)
labels[c(1,3,5,7,9)]<-c(paste0(round(length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$pLI>=0.9))/
                                       length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(universe_df$lethal_inheritance=="AR"&universe_df$pLI>=0.9))/
                                       length(which((universe_df$lethal_inheritance=="AR")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI>=0.9))/
                                       length(which((universe_df$lethal_inheritance=="AR,AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(universe_df$lethal_inheritance=="AD"&universe_df$pLI>=0.9))/
                                       length(which((universe_df$lethal_inheritance=="AD")&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"),
                        paste0(round(length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$pLI>=0.9))/
                                       length(which((grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance))&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))*100,0),"%"))
 
                               
l <- ggplot(dat=inh_lof_human_lethal_Bm, aes(x=variable, y=value, fill=mis_constraint))+
  geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse = TRUE))+
  scale_y_continuous(expand = c(0, 0),labels = percent) +bar_theme()+
  labs(y = "Percent constrained")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="none")+
  geom_text(aes(label = labels),size=3,position = "identity",vjust=-1)

ggarrange(m,l)
ggsave("output/Figures/2D.pdf",height=5, width=13, units='cm')
rm(e,inh_lof,inh_lof_human_lethal_B,inh_lof_human_lethal_Bm,inh_lofm,inh_mis,inh_mis_human_lethal_B,inh_mis_human_lethal_Bm,inh_mism,l,m,labelsf,labelsmis)

