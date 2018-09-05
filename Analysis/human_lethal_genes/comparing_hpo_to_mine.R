a<-which(grepl("Stillbirth",universe_df$hpo_names)|
      grepl("Neonatal death",universe_df$hpo_names)|
      grepl("Death in infancy",universe_df$hpo_names)|
      grepl("Death in childhood",universe_df$hpo_names)|
      grepl("Sudden death",universe_df$hpo_names)|
      grepl("Death in early adulthood",universe_df$hpo_names))

length(which(!is.na(universe_df$lethal_mouse[a])))

length(which(universe_df$human_lethal[a]=="Y"))
length(which(universe_df$human_lethal[a]=="N"))

hpo_lethal <- data.frame(gene=universe_df$gene[a],hpo_names=universe_df$hpo_names[a],my_list=universe_df$human_lethal[a])
lethal_terms <- c("Stillbirth","Neonatal death","Death in infancy","Death in childhood","Sudden death", "Death in early adulthood")
hpo_lethal$lethal_hit<- lapply(hpo_lethal$hpo_names,function(x) lethal_terms[which(lethal_terms%in%x)])
hpo_lethal<-hpo_lethal[,-2]
rm(a,lethal_terms)

#################graph of how many in each cat i caught#####################
counts<-data.frame(age.of.death=c("Stillbirth","Neonatal death","Death in infancy",
                                  "Death in childhood","Sudden death","Death in early adulthood"),
                   yes = c(length(which(grepl("Stillbirth",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="Y")),
                              length(which(grepl("Neonatal death",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="Y")),
                              length(which(grepl("Death in infancy",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="Y")),
                              length(which(grepl("Death in childhood",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="Y")),
                              length(which(grepl("Sudden death",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="Y")),
                              length(which(grepl("Death in early adulthood",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="Y"))),     
                   no= c(length(which(grepl("Stillbirth",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N")),
                   length(which(grepl("Neonatal death",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N")),
                   length(which(grepl("Death in infancy",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N")),
                   length(which(grepl("Death in childhood",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N")),
                   length(which(grepl("Death in early adulthood",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N")),
                   length(which(grepl("Sudden death",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N"))))
                              

counts$age.of.death<-factor(counts$age.of.death,levels=counts$age.of.death)
counts <- melt(counts)
names(counts)<-c("HPO_term","my_list","no_genes")

ggplot(data=counts,aes(x=HPO_term,y=no_genes,fill=my_list))+geom_bar(stat="identity")+coord_flip()
ggsave("hpo_comparison.pdf")



write.table(hpo_lethal$gene[which(hpo_lethal$lethal_hit=="Death in infancy"&hpo_lethal$my_list=="N")],"deathininfancy_notcaught.txt")
write.table(hpo_lethal$gene[which((grepl("Stillbirth",hpo_lethal$lethal_hit)|grepl("Neonatal death",hpo_lethal$lethal_hit))&hpo_lethal$my_list=="N")],"stillbirthneonatal_notcaught.txt")
write.table(hpo_lethal$gene[which(grepl("Death in childhood",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N")],"deathchildhood_notcaught.txt")
write.table(hpo_lethal$gene[which(grepl("Death in early adulthood",hpo_lethal$lethal_hit)&hpo_lethal$my_list=="N")],"deathadult_notcaught.txt")

