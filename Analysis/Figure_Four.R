source("plot_functions.R")

#4A:Scatter plot showing levels of genetic tolerance to LoF (pLI) or missense constraint ####
labelrec1<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI>=0.9&
                                       (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                          universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                     (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                        universe_df$lethal_inheritance=="MT,AR")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI>=0.9&
                                                                                                                                                                                      (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                                                                                                         universe_df$lethal_inheritance=="MT,AR"))))

labelrec2<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI>=0.9&
                                       (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                          universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                     (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                        universe_df$lethal_inheritance=="MT,AR")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI>=0.9&
                                                                                                                                                                                      (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                                                                                                         universe_df$lethal_inheritance=="MT,AR"))))
labelrec3<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI<0.9&
                                       (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                          universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                     (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                        universe_df$lethal_inheritance=="MT,AR")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI<0.9&
                                                                                                                                                                                      (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                                                                                                         universe_df$lethal_inheritance=="MT,AR"))))
labelrec4<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI<0.9&
                                       (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                          universe_df$lethal_inheritance=="MT,AR")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                     (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                        universe_df$lethal_inheritance=="MT,AR")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI<0.9&
                                                                                                                                                                                      (universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|
                                                                                                                                                                                         universe_df$lethal_inheritance=="MT,AR"))))

rec<-ggplot(universe_df[which(universe_df$omim=="Y"&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.05,alpha=0.2)+scatter_theme()+
  geom_point(data=universe_df[which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),], 
             aes(x=mis_z,y=pLI), colour="purple", size=0.2)+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("Recessive Human Lethal Genes \n (AR,XLr) \n n=",length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR"))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14)) +annotate(geom="text", x=-7, y=0.96,size=2, label=labelrec1)+annotate(geom="text", x=11, y=0.96,size=2, label=labelrec2)+annotate(geom="text", x=11, y=0.06,size=2, label=labelrec3)+annotate(geom="text", x=-7, y=0.06,size=2, label=labelrec4)


labeldom1<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI>=0.9&
                                    (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                       universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                  (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                     universe_df$lethal_inheritance=="MT,XLd")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI>=0.9&
                                                                                                                                                                                 (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                                                                                                   universe_df$lethal_inheritance=="MT,XLd"))))

labeldom2<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI>=0.9&
                                    (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                       universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                  (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                     universe_df$lethal_inheritance=="MT,XLd")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI>=0.9&
                                                                                                                                                                                (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                                                                                                   universe_df$lethal_inheritance=="MT,XLd"))))
labeldom3<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI<0.9&
                                    (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                       universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                  (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                     universe_df$lethal_inheritance=="MT,XLd")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09&universe_df$pLI<0.9&
                                                                                                                                                                                 (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                                                                                                    universe_df$lethal_inheritance=="MT,XLd"))))
labeldom4<-paste0(round(length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI<0.9&
                                    (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                       universe_df$lethal_inheritance=="MT,XLd")))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&
                                                                                                  (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                     universe_df$lethal_inheritance=="MT,XLd")))*100,0),"% \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09&universe_df$pLI<0.9&
                                                                                                                                                                                   (universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|
                                                                                                                                                                                      universe_df$lethal_inheritance=="MT,XLd"))))
dom<- ggplot(universe_df[which(universe_df$omim=="Y"&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,alpha=0.2)+scatter_theme()+
  geom_point(data=universe_df[which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),],
             aes(x=mis_z,y=pLI), colour="springgreen4", size=0.2)+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("Dominant Human Lethal Genes \n (AD,XLd) \n n=",length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd"))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14)) +annotate(geom="text", x=-7, y=0.96,size=2, label=labeldom1)+annotate(geom="text", x=11, y=0.96,size=2, label=labeldom2)+annotate(geom="text", x=11, y=0.06,size=2, label=labeldom3)+annotate(geom="text", x=-7, y=0.06,size=2, label=labeldom4)



ggarrange(rec,dom)
ggsave("output/Figures/4A.pdf",height=8.9, width=17.8, units='cm')
rm(rec,dom,labelrec1,labelrec2,labelrec3,labelrec4,labeldom1,labeldom2,labeldom3,labeldom4)




#4B:Table comparing constraint levels in OMIM non-lethal and OMIM lethal genes, separated by inheritance####
mt_no<-length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                       universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac)))
ar_no<-length(which((universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac))))
ar_ad_no<-length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac))))
ad_no<-length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac))))
xl_no<-length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac))))

mt_lethal_no<-length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&!is.na(universe_df$exac)))
ar_lethal_no<-length(which((universe_df$lethal_inheritance=="AR"&!is.na(universe_df$exac))))
ar_ad_lethal_no<-length(which((universe_df$lethal_inheritance=="AR,AD"&!is.na(universe_df$exac))))
ad_lethal_no<-length(which((universe_df$lethal_inheritance=="AD"&!is.na(universe_df$exac))))
xl_lethal_no<-length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                        universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&!is.na(universe_df$exac))))
#stats tests ####
OR_test <- function(a,b,c,d){
  matrix <- matrix(c(a,b,c,d),nrow=2,dimnames=list(metric1 = c("Y","N"),metric2 = c("Y","N")))
  test <- fisher.test(matrix,alternative = "two.sided")
  return(test)
}

# MT MISSENSE
mtmis <- OR_test(length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                        universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))+1,
        length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                        universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09))-1,
        length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z>=3.09))+1,
        length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z<3.09))-1)
# MT LOF
mtlof <- OR_test(length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                        universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))+1,
        length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                        universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))-1,
        length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$pLI>=0.9))+1,
        length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$pLI<0.9))-1)

# AR MISSENSE
armis <- OR_test(
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z<3.09)))
)

# AR LOF
arlof <- OR_test(
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI<0.9)))
)

# AR,AD MISSENSE
aradmis <- OR_test(
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z<3.09)))
)

# AR,AD LOF
aradlof <- OR_test(
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI<0.9)))  
)

# AD MISSENSE
admis <- OR_test(
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$mis_z<3.09))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09)))
)

# AD LOF
adlof <- OR_test(
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI<0.9))) 
)

# XL MISSENSE
xlmis <- OR_test(
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z>=3.09))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z<3.09)))
)

# XL LOF
xllof <- OR_test(
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI>=0.9))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI<0.9)))
)

#summary table ####
table <- data.frame(Inheritance = c("MT","AR","AR/AD","AD","XL"),
                    nonlethal_n = c(mt_no,ar_no,ar_ad_no,ad_no,xl_no),
                    perc_mis_intol = c(round(length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                                                              universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))/mt_no*100,0),
                                              round(length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/ar_no*100,0),
                                              round(length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/ar_ad_no*100,0),
                                              round(length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/ad_no*100,0),
                                              round(length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/xl_no*100,0)),
                    perc_lof_intol = c(round(length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                                                                     universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))/mt_no*100,0),
                                               round(length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/ar_no*100,0),
                                               round(length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/ar_ad_no*100,0),
                                               round(length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/ad_no*100,0),
                                               round(length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/xl_no*100,0)),
                    lethal_n = c(mt_lethal_no,ar_lethal_no,ar_ad_lethal_no,ad_lethal_no,xl_lethal_no),
                    perc_mis_intol = c(round(length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z>=3.09))/mt_lethal_no*100,0),
                                                    round(length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z>=3.09)))/ar_lethal_no*100,0),
                                                    round(length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z>=3.09)))/ar_ad_lethal_no*100,0),
                                                    round(length(which((universe_df$lethal_inheritance=="AD"&universe_df$mis_z>=3.09)))/ad_lethal_no*100,0),
                                                    round(length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                                                                     universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z>=3.09)))/xl_lethal_no*100,0)),
                    pval = c(round(mtmis$p.value,3),round(armis$p.value,3),round(aradmis$p.value,3),round(admis$p.value,3),round(xlmis$p.value,3)),
                    per_lof_intol = c(round(length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$pLI>=0.9))/mt_lethal_no*100,0),
                                               round(length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI>=0.9)))/ar_lethal_no*100,0),
                                               round(length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI>=0.9)))/ar_ad_lethal_no*100,0),
                                               round(length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI>=0.9)))/ad_lethal_no*100,0),
                                               round(length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                                                                      universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI>=0.9)))/xl_lethal_no*100,0)),
                    pval= c(round(mtlof$p.value,3),round(arlof$p.value,3),round(aradlof$p.value,3),round(adlof$p.value,3),round(xllof$p.value,3)) )

g <- tableGrob(table)
ggplot2::ggsave('output/Figures/4B.pdf',g,height=5, width=25, units='cm')

rm(table,g,mt_no,mt_lethal_no,ad_no,ad_lethal_no,ar_no,ar_lethal_no,ar_ad_lethal_no,ar_ad_no,armis,arlof,mtmis,mtlof,aradmis,aradlof,admis,adlof,xlmis,xllof,xl_no,xl_lethal_no)





#4C:OMIM prenatal/infantile lethality genes are strongly associated with lethal murine phenotypes####

mouse_lethal_prop <- data.frame(non_OMIM=c(length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))),
                                OMIM=c(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y")),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y")),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))),
                                Human_Lethal_B=c(length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y")),
                                                 length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y")),
                                                 length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))),
                                Human_Lethal_A=c(length(which(universe_df$human_lethal_A=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_A=="Y"&universe_df$mouse_ko=="Y")),
                                                 length(which(universe_df$human_lethal_A=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_A=="Y"&universe_df$mouse_ko=="Y")),
                                                 length(which(universe_df$human_lethal_A=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_A=="Y"&universe_df$mouse_ko=="Y"))))

levels<- c(paste0("non-OMIM \n n = ",length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))),
           paste0("OMIM \n n = ",length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))),
           paste0("Humanl Lethal \n List B \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))),
           paste0("Human Lethal \n List A \n n = ",length(which(universe_df$human_lethal_A=="Y"&universe_df$mouse_ko=="Y"))))

mouse_lethal_propm <- melt(mouse_lethal_prop)
mouse_lethal_propm$mis_constraint <- rep(c("Mouse Lethal     ","Premature death     ", "Mouse non-Lethal     "), 4)
mouse_lethal_propm$mis_constraint<-factor(mouse_lethal_propm$mis_constraint,levels=c("Mouse Lethal     ","Premature death     ", "Mouse non-Lethal     "))
mouse_lethal_propm$variable <- rep(levels,each=3)
mouse_lethal_propm$variable <- factor(mouse_lethal_propm$variable, levels = levels)

labels <- rep("",12)
labels[c(2,5,8,11)]<-c(paste0(round(length(which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|grepl("premature death",universe_df$all_MP_phen))))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))*100,0),"%"),
                       paste0(round(length(which(universe_df$omim=="Y"&(universe_df$lethal_mouse=="Y"|grepl("premature death",universe_df$all_MP_phen))))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))*100,0),"%"),
                       paste0(round(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_mouse=="Y"|grepl("premature death",universe_df$all_MP_phen))))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))*100,0),"%"),
                       paste0(round(length(which(universe_df$human_lethal_A=="Y"&(universe_df$lethal_mouse=="Y"|grepl("premature death",universe_df$all_MP_phen))))/length(which(universe_df$human_lethal_A=="Y"&universe_df$mouse_ko=="Y"))*100,0),"%"))


f <- ggplot(dat=mouse_lethal_propm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0),limits=c(0,1)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","grey","steelblue3"))+theme(legend.position="right")+geom_text(aes(label = labels),position = position_stack(reverse=TRUE,vjust=2))

ggsave("output/Figures/4C.pdf",height=10, width=17, units='cm')
rm(mouse_lethal_prop,mouse_lethal_propm,f,levels,labels)
