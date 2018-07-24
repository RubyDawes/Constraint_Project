
omim_ko <- data.frame(group=c("Available","Unavailable"),
                      value=c(length(which(universe_df$omim=="Y"&!is.na(universe_df$mouse_ko))),
                              length(which(universe_df$omim=="Y"&is.na(universe_df$mouse_ko)))))
nonomim_ko <- data.frame(group=c("Available","Unavailable"),
                      value=c(length(which(is.na(universe_df$omim)&!is.na(universe_df$mouse_ko))),
                              length(which(is.na(universe_df$omim)&is.na(universe_df$mouse_ko)))))
omim_ko$group <- factor(omim_ko$group, levels = omim_ko$group)
nonomim_ko$group <- factor(nonomim_ko$group, levels = nonomim_ko$group)
omim_ko_pie <- ggplot(omim_ko, aes(x="", y=value, fill=group))
omim_ko_pie<- omim_ko_pie+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
omim_ko_pie<- omim_ko_pie+scale_fill_manual(values=c("slategray1","slategray"))+theme(legend.position="none")
omim_ko_pie<- omim_ko_pie+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))
omim_ko_pie<- omim_ko_pie+ggtitle(paste("OMIM genes \n n=",length(which(universe_df$omim=="Y"))))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))

nonomim_ko_pie <- ggplot(nonomim_ko, aes(x="", y=value, fill=group))
nonomim_ko_pie<- nonomim_ko_pie+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
nonomim_ko_pie<- nonomim_ko_pie+scale_fill_manual(values=c("slategray1","slategray"))+theme(legend.position="none")
nonomim_ko_pie<- nonomim_ko_pie+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))
nonomim_ko_pie<- nonomim_ko_pie+ggtitle(paste("Non-OMIM genes \n n=",length(which(is.na(universe_df$omim)))))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))


ggarrange(omim_ko_pie,nonomim_ko_pie,common.legend=TRUE,legend="bottom")
  ggsave("Sandra_Figures/Figs/fig2a-KOpies.pdf",height=6, width=17.8, units='cm')
  dev.off()
rm(omim_ko,nonomim_ko,omim_ko_pie,nonomim_ko_pie)