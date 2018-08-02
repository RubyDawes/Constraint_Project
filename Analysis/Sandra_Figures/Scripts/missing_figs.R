#mouse lethal genes more likely to be oMIM genes
leth_omim <- data.frame(non_OMIM=c(length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&!is.na(universe_df$lethal_mouse))),
                               length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"))/length(which(is.na(universe_df$omim)&!is.na(universe_df$lethal_mouse)))),
                      OMIM=c(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse))),
                             length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)))))
                      
leth_omimm <- melt(leth_omim)
leth_omimm$mouse_leth <- rep(c("Lethal     ", "Non-Lethal     "), 2)

f <- ggplot(dat=leth_omimm, aes(x=variable, y=value, fill=mouse_leth))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig-mouseleth_prop.pdf",height=7, width=7, units='cm')
rm(leth_omim,leth_omimm,f)

#cell essential genes more likely to be oMIM genes
celless_omim <- data.frame(non_OMIM=c(length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y"))/length(which(is.na(universe_df$omim)&!is.na(universe_df$cell_essential))),
                                   length(which(is.na(universe_df$omim)&universe_df$cell_essential=="N"))/length(which(is.na(universe_df$omim)&!is.na(universe_df$cell_essential)))),
                        OMIM=c(length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$cell_essential))),
                               length(which(universe_df$omim=="Y"&universe_df$cell_essential=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$cell_essential)))))

celless_omimm <- melt(celless_omim)
celless_omimm$mouse_leth <- rep(c("Essential     ", "Non-Essential     "), 2)

f <- ggplot(dat=celless_omimm, aes(x=variable, y=value, fill=mouse_leth))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("sandybrown","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/fig-celless_prop.pdf",height=7, width=7, units='cm')
rm(celless_omim,celless_omimm,f)
