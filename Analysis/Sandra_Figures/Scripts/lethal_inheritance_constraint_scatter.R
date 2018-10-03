
ggplot(universe_df[which(universe_df$omim=="Y"&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,alpha=0.3)+scatter_theme_goodsizes()+
  geom_point(data=universe_df[which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),], aes(x=mis_z,y=pLI), colour="blue", size=1.5,alpha=0.7)+
  geom_point(data=universe_df[which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),], aes(x=mis_z,y=pLI), colour="red", size=1.5,alpha=0.7)+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("OMIM genes \n n=",length(which(universe_df$omim=="Y"&!is.na(universe_df$exac))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14)) 


rec<-ggplot(universe_df[which(universe_df$omim=="Y"&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.05,alpha=0.2)+scatter_theme_goodsizes()+
  geom_point(data=universe_df[which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),], 
             aes(x=mis_z,y=pLI), colour="purple", size=0.2)+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("Recessive Human Lethal Genes \n (AR,XLr) \n n=",length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR"))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14)) 


dom<- ggplot(universe_df[which(universe_df$omim=="Y"&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,alpha=0.2)+scatter_theme_goodsizes()+
  geom_point(data=universe_df[which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),],
             aes(x=mis_z,y=pLI), colour="springgreen4", size=0.2)+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("Dominant Human Lethal Genes \n (AD,XLd) \n n=",length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd"))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14)) 


ggarrange(rec,dom)
ggsave("Analysis/Sandra_Figures/Figs/lethal_inh_constraint_scatter.pdf",height=6, width=10, units='cm')



length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI<0.9&
               (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                  universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                             (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                universe_df$lethal_inheritance=="MT,XLd")))*100

length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI>=0.9&
               (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                  universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                             (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                universe_df$lethal_inheritance=="MT,XLd")))*100


length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI>=0.9&
               (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                  universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                             (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                universe_df$lethal_inheritance=="MT,XLd")))*100

length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI<0.9&
               (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                  universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                             (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                universe_df$lethal_inheritance=="MT,XLd")))*100



length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI<0.9&
               (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                  universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                             (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                universe_df$lethal_inheritance=="MT,AR")))*100

length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI>=0.9&
               (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                  universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                            (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                               universe_df$lethal_inheritance=="MT,AR")))*100

length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI>=0.9&
               (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                  universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                            (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                               universe_df$lethal_inheritance=="MT,AR")))*100

length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI<0.9&
               (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                  universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                            (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                               universe_df$lethal_inheritance=="MT,AR")))*100
