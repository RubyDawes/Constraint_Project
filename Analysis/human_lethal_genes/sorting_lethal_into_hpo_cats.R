#adding column for age of death
load("output/Data/human_lethal_genes.rda")

#restructuring matches column for easier searching
for (i in which(lengths(lethal_genes$matches)>1)) {
  lethal_genes$matches[[i]] <- paste(lethal_genes$matches[[i]][1:length(lethal_genes$matches[[i]])],collapse=',')
}
lethal_genes$matches<-unlist(lapply(lethal_genes$matches,function(x) gsub('\"', "", x, fixed = TRUE)))
##########sorting lethal genesinto age of death#########
lethal_genes$prenatal <- ifelse((  grepl('death in utero~2',lethal_genes$matches)|
                                     grepl('delivered dead~3',lethal_genes$matches)|
                                     grepl('lethal before birth~4',lethal_genes$matches)|
                                     grepl('lethal in utero~2',lethal_genes$matches)|
                                     grepl('die in utero~2',lethal_genes$matches)|
                                     grepl('spontaneous abortion',lethal_genes$matches)|
                                     grepl('fetal demise~3',lethal_genes$matches)|
                                     grepl('embryonic lethal~3',lethal_genes$matches)|
                                     grepl('embryonically lethal~3',lethal_genes$matches)|
                                     grepl('prenatally lethal~3',lethal_genes$matches)|
                                     grepl('fetal death',lethal_genes$matches)|
                                     grepl('stillborn fetus~3',lethal_genes$matches)|
                                     grepl('lethal before birth~4',lethal_genes$matches)|
                                     grepl('lethal prenatal~4',lethal_genes$matches)),"Y","N")
lethal_genes$neonatal <- ifelse((  grepl('death occurred hours~4',lethal_genes$matches)|
                                     grepl('died days~3',lethal_genes$matches)|
                                     grepl('died shortly after birth~2',lethal_genes$matches)|
                                     grepl('neonatally lethal~2',lethal_genes$matches)|
                                     grepl('lethal neonatal~3',lethal_genes$matches)|
                                     grepl('died perinatal~4',lethal_genes$matches)|
                                     grepl('born died shortly after~2',lethal_genes$matches)|
                                     grepl('death neonatal~4',lethal_genes$matches)|
                                     grepl('dead perinatal~4',lethal_genes$matches)|
                                     grepl('death perinatal~4',lethal_genes$matches)|
                                     grepl('dead perinatal~4',lethal_genes$matches)|
                                     grepl('lethal perinatal~4',lethal_genes$matches)|
                                     grepl('lethal neonatal~3',lethal_genes$matches)|
                                     grepl('died hours after delivery~5',lethal_genes$matches)|
                                     grepl('died hours after birth~7',lethal_genes$matches)|
                                     grepl('died days~4',lethal_genes$matches)|
                                     grepl('death at days~3',lethal_genes$matches)|
                                     grepl('died at hours~2',lethal_genes$matches)|
                                     grepl('died days of life~3',lethal_genes$matches)|
                                     grepl('died hours after birth~5',lethal_genes$matches)),"Y","N")
lethal_genes$infancy <- ifelse((  grepl('death first year of life~3',lethal_genes$matches)|
  grepl('death months life~4',lethal_genes$matches)|
  grepl('died months life~4',lethal_genes$matches)|
  grepl('dead months life~4',lethal_genes$matches)|
  grepl('lethal first year of life~3',lethal_genes$matches)|
  grepl('Died at weeks~3',lethal_genes$matches)|
  grepl('died at months~4',lethal_genes$matches)|
  grepl('infants die~3',lethal_genes$matches)|
  grepl('die within the first year of life',lethal_genes$matches)|
  grepl('the infant died',lethal_genes$matches)|
  grepl('died infancy~3',lethal_genes$matches)|
  grepl('died infants~3',lethal_genes$matches)|
  grepl('died as infants',lethal_genes$matches)|
  grepl('death infancy~4',lethal_genes$matches)|
  grepl('lethal infantile~4',lethal_genes$matches)|
  grepl('lethal first year of life~3',lethal_genes$matches)|
  grepl('infants died~3',lethal_genes$matches)|
  grepl('died 1 year~3',lethal_genes$matches)|
  grepl('died at months of age~3',lethal_genes$matches)|
  grepl('died in infancy~4',lethal_genes$matches)|   
  grepl('death infant~4',lethal_genes$matches)|
  grepl('died at months~3',lethal_genes$matches)|
  grepl('death age months~4',lethal_genes$matches)|
  grepl('death by months~5',lethal_genes$matches)|
  grepl('died within months~2',lethal_genes$matches)|
  grepl('infant fatal~3',lethal_genes$matches)|   
  grepl('infant died',lethal_genes$matches)|
  grepl('fatal infancy~3',lethal_genes$matches)),"Y","N")

lethal_genes$sort <- rep(NA,length(lethal_genes$gene))
lethal_genes$sort[which(lethal_genes$prenatal=="Y")] <- "prenatal"
lethal_genes$sort[which(lethal_genes$neonatal == "Y")]<- ifelse(is.na(lethal_genes$sort[which(lethal_genes$neonatal == "Y")]),"neonatal",paste(lethal_genes$sort[which(lethal_genes$neonatal == "Y")],"neonatal",sep=","))
lethal_genes$sort[which(lethal_genes$infancy == "Y")]<- ifelse(is.na(lethal_genes$sort[which(lethal_genes$infancy == "Y")]),"infancy",paste(lethal_genes$sort[which(lethal_genes$infancy == "Y")],"infancy",sep=","))
############ venn of overlap of neonatal, prenatal, infancy genes #############
overrideTriple<- "ye"
png('Analysis/Sandra_Figures/Figs/venn_sorted_lethal_genes.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(lethal_genes$prenatal=="Y")), area2 = length(which(lethal_genes$neonatal=="Y")), area3 = length(which(lethal_genes$infancy=="Y")), 
                 n12 = length(which(lethal_genes$prenatal=="Y"&lethal_genes$neonatal=="Y")),
                 n23 = length(which(lethal_genes$infancy=="Y"&lethal_genes$neonatal=="Y")), 
                 n13 = length(which(lethal_genes$prenatal=="Y"&lethal_genes$infancy=="Y")),
                 n123 = length(which(lethal_genes$prenatal=="Y"&lethal_genes$neonatal=="Y"&lethal_genes$infancy=="Y")), category = c("Prenatal", "Neonatal", "Infancy"), lty = "blank", 
                 fill = c("#f1c266", "#f166c5", "#66f1c0"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(4,3))
dev.off()

############how many mouse lethal in each category###############
death_phens_mouse <- data.frame("prenatal"=c(length(which(lethal_genes$prenatal=="Y"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$prenatal=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$prenatal=="Y"&lethal_genes$lethal_mouse=="N"))/length(which(lethal_genes$prenatal=="Y"&!is.na(lethal_genes$lethal_mouse)))),
                                "neonatal"=c(length(which(lethal_genes$neonatal=="Y"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$neonatal=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$neonatal=="Y"&lethal_genes$lethal_mouse=="N"))/length(which(lethal_genes$neonatal=="Y"&!is.na(lethal_genes$lethal_mouse)))),
                                "infancy"=c(length(which(lethal_genes$infancy=="Y"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$infancy=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                            length(which(lethal_genes$infancy=="Y"&lethal_genes$lethal_mouse=="N"))/length(which(lethal_genes$infancy=="Y"&!is.na(lethal_genes$lethal_mouse)))),
                                "unsorted"=c(length(which(is.na(lethal_genes$sort)&lethal_genes$lethal_mouse=="Y"))/length(which(is.na(lethal_genes$sort)&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(is.na(lethal_genes$sort)&lethal_genes$lethal_mouse=="N"))/length(which(is.na(lethal_genes$sort)&!is.na(lethal_genes$lethal_mouse)))))
      
death_phens_mousem <- melt(death_phens_mouse)
death_phens_mousem$lethality <- rep(c("Lethal in a mouse", "Non-lethal in a mouse"), 4)

d <- ggplot(dat=death_phens_mousem, aes(x=variable, y=value, fill=lethality))
d<- d+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
d<- d+labs(x = "HPO term")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="bottom")+coord_flip()

png("output/figures/sorted_lethal_proportions_phens.png",width=1500,height=1000,type="quartz",res=150,bg = "transparent")
d
dev.off()


############how many mouse lethal in each category including prem death in grey###############
death_phens_mouse <- data.frame("prenatal"=c(length(which(lethal_genes$prenatal=="Y"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$prenatal=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$prenatal=="Y"&lethal_genes$lethal_mouse=="N"&grepl("premature death",lethal_genes$all_MP_phen)))/length(which(lethal_genes$prenatal=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$prenatal=="Y"&lethal_genes$lethal_mouse=="N"&!grepl("premature death",lethal_genes$all_MP_phen)))/length(which(lethal_genes$prenatal=="Y"&!is.na(lethal_genes$lethal_mouse)))),
                                "neonatal"=c(length(which(lethal_genes$neonatal=="Y"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$neonatal=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$neonatal=="Y"&lethal_genes$lethal_mouse=="N"&grepl("premature death",lethal_genes$all_MP_phen)))/length(which(lethal_genes$neonatal=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$neonatal=="Y"&lethal_genes$lethal_mouse=="N"&!grepl("premature death",lethal_genes$all_MP_phen)))/length(which(lethal_genes$neonatal=="Y"&!is.na(lethal_genes$lethal_mouse)))),
                                "infancy"=c(length(which(lethal_genes$infancy=="Y"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$infancy=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                            length(which(lethal_genes$infancy=="Y"&lethal_genes$lethal_mouse=="N"&grepl("premature death",lethal_genes$all_MP_phen)))/length(which(lethal_genes$infancy=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                            length(which(lethal_genes$infancy=="Y"&lethal_genes$lethal_mouse=="N"&!grepl("premature death",lethal_genes$all_MP_phen)))/length(which(lethal_genes$infancy=="Y"&!is.na(lethal_genes$lethal_mouse)))),
                                "unsorted"=c(length(which(is.na(lethal_genes$sort)&lethal_genes$lethal_mouse=="Y"))/length(which(is.na(lethal_genes$sort)&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(is.na(lethal_genes$sort)&lethal_genes$lethal_mouse=="N"&grepl("premature death",lethal_genes$all_MP_phen)))/length(which(is.na(lethal_genes$sort)&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(is.na(lethal_genes$sort)&lethal_genes$lethal_mouse=="N"&!grepl("premature death",lethal_genes$all_MP_phen)))/length(which(is.na(lethal_genes$sort)&!is.na(lethal_genes$lethal_mouse)))))

      
death_phens_mousem <- melt(death_phens_mouse)
death_phens_mousem$lethality <- rep(c("a-Lethal in a mouse", "b-premature death in mouse","c-Non-lethal in a mouse"), 4)

d <- ggplot(dat=death_phens_mousem, aes(x=variable, y=value, fill=lethality))
d<- d+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
d<- d+labs(x = "HPO term")+scale_fill_manual(values=c("black","grey","steelblue3"))+theme(legend.position="bottom")+coord_flip()
      
      png("output/figures/sorted_lethal_proportions_phens_withpremdeath.png",width=1500,height=1000,type="quartz",res=150,bg = "transparent")
      d
      dev.off()
      
      
###########how many mouse lethal but only taking earliest death category############
death_phens_mouse <- data.frame("prenatal"=c(length(which(lethal_genes$prenatal=="Y"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$prenatal=="Y"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$prenatal=="Y"&lethal_genes$lethal_mouse=="N"))/length(which(lethal_genes$prenatal=="Y"&!is.na(lethal_genes$lethal_mouse)))),
                                "neonatal"=c(length(which(lethal_genes$neonatal=="Y"&lethal_genes$prenatal=="N"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$neonatal=="Y"&lethal_genes$prenatal=="N"&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(lethal_genes$neonatal=="Y"&lethal_genes$prenatal=="N"&lethal_genes$lethal_mouse=="N"))/length(which(lethal_genes$neonatal=="Y"&lethal_genes$prenatal=="N"&!is.na(lethal_genes$lethal_mouse)))),
                                "infancy"=c(length(which(lethal_genes$infancy=="Y"&lethal_genes$prenatal=="N"&lethal_genes$neonatal=="N"&lethal_genes$lethal_mouse=="Y"))/length(which(lethal_genes$infancy=="Y"&lethal_genes$prenatal=="N"&lethal_genes$neonatal=="N"&!is.na(lethal_genes$lethal_mouse))),
                                            length(which(lethal_genes$infancy=="Y"&lethal_genes$prenatal=="N"&lethal_genes$neonatal=="N"&lethal_genes$lethal_mouse=="N"))/length(which(lethal_genes$infancy=="Y"&lethal_genes$prenatal=="N"&lethal_genes$neonatal=="N"&!is.na(lethal_genes$lethal_mouse)))),
                                "unsorted"=c(length(which(is.na(lethal_genes$sort)&lethal_genes$lethal_mouse=="Y"))/length(which(is.na(lethal_genes$sort)&!is.na(lethal_genes$lethal_mouse))),
                                             length(which(is.na(lethal_genes$sort)&lethal_genes$lethal_mouse=="N"))/length(which(is.na(lethal_genes$sort)&!is.na(lethal_genes$lethal_mouse)))))


death_phens_mousem <- melt(death_phens_mouse)
death_phens_mousem$lethality <- rep(c("Lethal in a mouse", "Non-lethal in a mouse"), 4)

d <- ggplot(dat=death_phens_mousem, aes(x=variable, y=value, fill=lethality))
d<- d+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
d<- d+labs(x = "HPO term")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="bottom")+coord_flip()

png("output/figures/sorted_lethal_proportions_phens_onlyearliest.png",width=1500,height=1000,type="quartz",res=150,bg = "transparent")
d
dev.off()
#no real difference- slight decrease in infancy but eh

len



#in the process
##############not used###########
# died before the age death within birth~4    death after birth~3 
#lethal early life~4 lethal malfmation~3  no development milestones~2 severe lethal~2  early lethal~2  die after birth~2  lethal disorder~4  
# lethal disorder~4  neonatal severe~3 death early life~3  early death died early age~3 


##############Death before birth############
lethal_genes$matches[which(
  grepl('death in utero~2',lethal_genes$matches)|
    grepl('delivered dead~3',lethal_genes$matches)|
    grepl('lethal before birth~4',lethal_genes$matches)|
    grepl('lethal in utero~2',lethal_genes$matches)|
    grepl('die in utero~2',lethal_genes$matches)|
    grepl('spontaneous abortion',lethal_genes$matches)|
    grepl('fetal demise~3',lethal_genes$matches)|
    grepl('embryonic lethal~3',lethal_genes$matches)|
    grepl('embryonically lethal~3',lethal_genes$matches)|
    grepl('prenatally lethal~3',lethal_genes$matches)|
    grepl('fetal death',lethal_genes$matches)|
    grepl('stillborn fetus~3',lethal_genes$matches)|
    grepl('lethal before birth~4',lethal_genes$matches)|
    grepl('lethal prenatal~4',lethal_genes$matches)
)]  



#Death days after birth- perinatal/neonatal
lethal_genes$matches[which(
  grepl('death occurred hours~4',lethal_genes$matches)|
    grepl('died days~3',lethal_genes$matches)|
    grepl('died shortly after birth~2',lethal_genes$matches)|
    grepl('neonatally lethal~2',lethal_genes$matches)|
    grepl('lethal neonatal~3',lethal_genes$matches)|
    grepl('died perinatal~4',lethal_genes$matches)|
    grepl('born died shortly after~2',lethal_genes$matches)|
    grepl('death neonatal~4',lethal_genes$matches)|
    grepl('dead perinatal~4',lethal_genes$matches)|
    grepl('death perinatal~4',lethal_genes$matches)|
    grepl('dead perinatal~4',lethal_genes$matches)|
    grepl('lethal perinatal~4',lethal_genes$matches)|
    grepl('lethal neonatal~3',lethal_genes$matches)|
    grepl('died hours after delivery~5',lethal_genes$matches)|
    grepl('died hours after birth~7',lethal_genes$matches)|
    grepl('died days~4',lethal_genes$matches)|
    grepl('death at days~3',lethal_genes$matches)|
    grepl('died at hours~2',lethal_genes$matches)|
    grepl('died days of life~3',lethal_genes$matches)|
    grepl('died hours after birth~5',lethal_genes$matches)
)]    


 


 
#death in infancy
lethal_genes$matches[which(
  grepl('death first year of life~3',lethal_genes$matches)|
    grepl('death months life~4',lethal_genes$matches)|
    grepl('died months life~4',lethal_genes$matches)|
    grepl('dead months life~4',lethal_genes$matches)|
    grepl('lethal first year of life~3',lethal_genes$matches)|
    grepl('Died at weeks~3',lethal_genes$matches)|
    grepl('died at months~4',lethal_genes$matches)|
    grepl('infants die~3',lethal_genes$matches)|
    grepl('die within the first year of life',lethal_genes$matches)|
    grepl('the infant died',lethal_genes$matches)|
    grepl('died infancy~3',lethal_genes$matches)|
    grepl('died infants~3',lethal_genes$matches)|
    grepl('died as infants',lethal_genes$matches)|
    grepl('death infancy~4',lethal_genes$matches)|
    grepl('lethal infantile~4',lethal_genes$matches)|
    grepl('lethal first year of life~3',lethal_genes$matches)|
    grepl('infants died~3',lethal_genes$matches)|
    grepl('died 1 year~3',lethal_genes$matches)|
    grepl('died at months of age~3',lethal_genes$matches)|
    grepl('died in infancy~4',lethal_genes$matches)|   
    grepl('death infant~4',lethal_genes$matches)|
    grepl('died at months~3',lethal_genes$matches)|
    grepl('death age months~4',lethal_genes$matches)|
    grepl('death by months~5',lethal_genes$matches)|
    grepl('died within months~2',lethal_genes$matches)|
    grepl('infant fatal~3',lethal_genes$matches)|   
    grepl('infant died',lethal_genes$matches)|
    grepl('fatal infancy~3',lethal_genes$matches)
)]


 
