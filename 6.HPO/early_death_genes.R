
death_phens_mouse <- data.frame("Stillbirth"=c(length(which(grepl("Stillbirth",universe_df$hpo_names)&universe_df$lethal_mouse=="Y"))/length(which(grepl("Stillbirth",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse))),
                                               length(which(grepl("Stillbirth",universe_df$hpo_names)&universe_df$lethal_mouse=="N"))/length(which(grepl("Stillbirth",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse)))),
                                "Neonatal death"=c(length(which(grepl("Neonatal death",universe_df$hpo_names)&universe_df$lethal_mouse=="Y"))/length(which(grepl("Neonatal death",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse))),
                                               length(which(grepl("Neonatal death",universe_df$hpo_names)&universe_df$lethal_mouse=="N"))/length(which(grepl("Neonatal death",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse)))),
                                "Death in infancy"=c(length(which(grepl("Death in infancy",universe_df$hpo_names)&universe_df$lethal_mouse=="Y"))/length(which(grepl("Death in infancy",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse))),
                                               length(which(grepl("Death in infancy",universe_df$hpo_names)&universe_df$lethal_mouse=="N"))/length(which(grepl("Death in infancy",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse)))),
                                "Death in childhood"=c(length(which(grepl("Death in childhood",universe_df$hpo_names)&universe_df$lethal_mouse=="Y"))/length(which(grepl("Death in childhood",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse))),
                                               length(which(grepl("Death in childhood",universe_df$hpo_names)&universe_df$lethal_mouse=="N"))/length(which(grepl("Death in childhood",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse)))),
                                "Sudden death"=c(length(which(grepl("Sudden death",universe_df$hpo_names)&universe_df$lethal_mouse=="Y"))/length(which(grepl("Sudden death",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse))),
                                               length(which(grepl("Sudden death",universe_df$hpo_names)&universe_df$lethal_mouse=="N"))/length(which(grepl("Sudden death",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse)))),
                                "Death in early adulthood"=c(length(which(grepl("Death in early adulthood",universe_df$hpo_names)&universe_df$lethal_mouse=="Y"))/length(which(grepl("Death in early adulthood",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse))),
                                               length(which(grepl("Death in early adulthood",universe_df$hpo_names)&universe_df$lethal_mouse=="N"))/length(which(grepl("Death in early adulthood",universe_df$hpo_names)&!is.na(universe_df$lethal_mouse)))))
                                
                                
death_phens_mousem <- melt(death_phens_mouse)
death_phens_mousem$lethality <- rep(c("Lethal in a mouse", "Non-lethal in a muose"), 6)

d <- ggplot(dat=death_phens_mousem, aes(x=variable, y=value, fill=lethality))
d<- d+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
d<- d+labs(x = "HPO term")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="bottom")+coord_flip()
d<- d+scale_y_continuous(breaks = pretty(death_phens_mousem$value, n = 5))
png("output/figures/lethal_proportions_phens.png",width=1500,height=1000,type="quartz",res=150,bg = "transparent")
d
dev.off()

length(which(grepl("Stillbirth",universe_df$hpo_names)&universe_df$constrained=="Y"))
length(which(grepl("Stillbirth",universe_df$hpo_names)&universe_df$constrained=="N"))
length(which(grepl("Neonatal death",universe_df$hpo_names)&universe_df$constrained=="Y"))
length(which(grepl("Neonatal death",universe_df$hpo_names)&universe_df$constrained=="N"))
length(which(grepl("Death in infancy",universe_df$hpo_names)&universe_df$constrained=="Y"))
length(which(grepl("Death in infancy",universe_df$hpo_names)&universe_df$constrained=="N"))
length(which(grepl("Death in childhood",universe_df$hpo_names)&universe_df$constrained=="Y"))
length(which(grepl("Death in childhood",universe_df$hpo_names)&universe_df$constrained=="N"))
length(which(grepl("Sudden death",universe_df$hpo_names)&universe_df$constrained=="Y"))
length(which(grepl("Sudden death",universe_df$hpo_names)&universe_df$constrained=="N"))
length(which(grepl("Death in early adulthood",universe_df$hpo_names)&universe_df$constrained=="Y"))
length(which(grepl("Death in early adulthood",universe_df$hpo_names)&universe_df$constrained=="N"))

8/(8+17)
6/(6+11)
32/(32+100)
5/(5+17)
5/(5+14)
9/(9+20)


