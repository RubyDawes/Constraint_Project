inh_leth <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="Y"))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$lethal_mouse))),
                           length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="N"))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$lethal_mouse)))),
                      AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$lethal_mouse))),
                           length(which(universe_df$Inheritance_pattern=="AR"&universe_df$lethal_mouse=="N"))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$lethal_mouse)))),
                      "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$lethal_mouse))),
                               length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$lethal_mouse=="N"))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$lethal_mouse)))),
                      AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$lethal_mouse))),
                           length(which(universe_df$Inheritance_pattern=="AD"&universe_df$lethal_mouse=="N"))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$lethal_mouse)))),
                      XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="Y"))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$lethal_mouse))),
                           length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="N"))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$lethal_mouse)))))

inh_lethm <- melt(inh_leth)
inh_lethm$mis_constraint <- rep(c("Lethal in a mouse     ", "Non-lethal in a mouse     "), 5)

f <- ggplot(dat=inh_lethm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="bottom")
ggsave("Sandra_Figures/Figs/lethality_vs_inheritance.pdf",height=7, width=7, units='cm')
