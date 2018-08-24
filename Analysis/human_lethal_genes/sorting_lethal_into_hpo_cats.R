######################adding column for age of death#####################
load("output/Data/human_lethal_genes.rda")

#restructuring matches column for easier searching
for (i in which(lengths(lethal_genes$matches)>1)) {
  lethal_genes$matches[[i]] <- paste(lethal_genes$matches[[i]][1:length(lethal_genes$matches[[i]])],collapse=',')
}
lethal_genes$matches<-unlist(lapply(lethal_genes$matches,function(x) gsub('\"', "", x, fixed = TRUE)))





##############not used###########
#the infant died fatal infancy~3 died befe the age died infancy~3 died infants~3death within birth~4 died as infants  bn died shtly after~2  death infancy~4 death infant~4 death after birth~3 
#death neonatal~4 died perinatal~4 dead perinatal~4 death perinatal~4 died perinatal~4 dead perinatal~4 
#lethal infantile~4 lethal perinatal~4 lethal prenatal~4 lethal early life~4 lethal malfmation~3lethal neonatal~3 no development milestones~2 severe lethal~2  early lethal~2  die after birth~2  lethal disorder~4  

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
    grepl('prenatally lethal~3',lethal_genes$matches)
)]  



#Death days after birth
lethal_genes$matches[which(
  grepl('death occurred hours~4',lethal_genes$matches)|
    grepl('died days~3',lethal_genes$matches)|
    grepl('neonatally lethal~2',lethal_genes$matches)|
    grepl('died shortly after birth~2',lethal_genes$matches)|
    grepl('died hours after birth~5',lethal_genes$matches)
)]    


#death weeks/months after birth
lethal_genes$matches[which(
  grepl('death first year of life~3',lethal_genes$matches)|
    grepl('death months life~4',lethal_genes$matches)|
    grepl('died months life~4',lethal_genes$matches)|
    grepl('dead months life~4',lethal_genes$matches)|
    grepl('lethal first year of life~3',lethal_genes$matches)|
    grepl('Died at weeks~3',lethal_genes$matches)|
    grepl('died at months~4',lethal_genes$matches)|
    grepl('infants die~3',lethal_genes$matches)|
    grepl('die within the first year of life',lethal_genes$matches)
)]



#death years after birth
lethal_genes$matches[which(
  grepl('die early childhood~3',lethal_genes$matches)
  
)] 

