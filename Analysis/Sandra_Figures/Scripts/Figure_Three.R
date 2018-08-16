#3a: inheritance stacked bar charts NMD and human lethal
inh_mis_nmd <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$mis_z>=3.09))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$mis_z<3.09))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                                   length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$mis_z>=3.09))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$mis_z<3.09))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac)))))

inh_mis_nmdm <- melt(inh_mis_nmd)
inh_mis_nmdm$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_mis_nmdm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig3ai.pdf",height=7, width=7, units='cm')
rm(inh_mis_nmd,inh_mis_nmdm,f)
inh_mis_human_lethal <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$mis_z>=3.09))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$mis_z<3.09))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                            length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$mis_z>=3.09))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$mis_z<3.09))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))))
inh_mis_human_lethalm <- melt(inh_mis_human_lethal)
inh_mis_human_lethalm$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_mis_human_lethalm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig3aii.pdf",height=7, width=7, units='cm')
rm(inh_mis_human_lethal,inh_mis_human_lethalm,f)

inh_lof_nmd <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$pLI>=0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$pLI<0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                                   length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$nmd=="Y"&!is.na(universe_df$exac)))),
                          XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$pLI>=0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac))),
                               length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&universe_df$pLI<0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$nmd=="Y"&!is.na(universe_df$exac)))))

inh_lof_nmdm <- melt(inh_lof_nmd)
inh_lof_nmdm$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_lof_nmdm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig3aiii.pdf",height=7, width=7, units='cm')
rm(inh_lof_nmd,inh_lof_nmdm,f)

inh_lof_human_lethal <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$pLI>=0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$pLI<0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                            length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))),
                                   XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$pLI>=0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&universe_df$pLI<0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal=="Y"&!is.na(universe_df$exac)))))

inh_lof_human_lethalm <- melt(inh_lof_human_lethal)
inh_lof_human_lethalm$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_lof_human_lethalm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig3aiv.pdf",height=7, width=7, units='cm')
rm(inh_lof_human_lethal,inh_lof_human_lethalm,f)

#Figure 3B proportions of mouse lethal for OMIM, NMD and human_lethal
mouse_lethal_prop <- data.frame(OMIM=c(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse))),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)))),
                                Human_Lethal=c(length(which(universe_df$human_lethal=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal=="Y"&!is.na(universe_df$lethal_mouse))),
                                               length(which(universe_df$human_lethal=="Y"&universe_df$lethal_mouse=="N"))/length(which(universe_df$human_lethal=="Y"&!is.na(universe_df$lethal_mouse)))))


mouse_lethal_propm <- melt(mouse_lethal_prop)
mouse_lethal_propm$mis_constraint <- rep(c("Mouse Lethal     ", "Mouse non-Lethal     "), 2)

f <- ggplot(dat=mouse_lethal_propm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig3bi.pdf",height=7, width=7, units='cm')
rm(mouse_lethal_prop,mouse_lethal_propm,f)

#adding in non-OMIM
mouse_lethal_prop <- data.frame(non_OMIM=c(length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&!is.na(universe_df$lethal_mouse))),
                                           length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"))/length(which(is.na(universe_df$omim)&!is.na(universe_df$lethal_mouse)))),
                                OMIM=c(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse))),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)))),
                                Human_Lethal=c(length(which(universe_df$human_lethal=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal=="Y"&!is.na(universe_df$lethal_mouse))),
                                               length(which(universe_df$human_lethal=="Y"&universe_df$lethal_mouse=="N"))/length(which(universe_df$human_lethal=="Y"&!is.na(universe_df$lethal_mouse)))))

mouse_lethal_propm <- melt(mouse_lethal_prop)
mouse_lethal_propm$mis_constraint <- rep(c("Mouse Lethal     ", "Mouse non-Lethal     "), 3)

f <- ggplot(dat=mouse_lethal_propm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig3bi-withnonomim.pdf",height=7, width=7, units='cm')
rm(mouse_lethal_prop,mouse_lethal_propm,f)


#Figure 3C HPO lethal annotations human lethal- mouse lethality

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
ggsave("Analysis/Sandra_Figures/Figs/fig3bii.pdf",height=7, width=14, units='cm')
