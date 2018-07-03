
omim_ko <- data.frame(group=c("Mouse Phenotyping data available","Mouse Phenotyping data unavailable"),
                      value=c(length(which(universe_df$omim=="Y"&!is.na(universe_df$mouse_ko))),
                              length(which(universe_df$omim=="Y"&is.na(universe_df$mouse_ko)))))
nonomim_ko <- data.frame(group=c("Mouse Phenotyping data available","Mouse Phenotyping data unavailable"),
                      value=c(length(which(is.na(universe_df$omim)&!is.na(universe_df$mouse_ko))),
                              length(which(is.na(universe_df$omim)&is.na(universe_df$mouse_ko)))))
omim_ko$group <- factor(omim_ko$group, levels = omim_ko$group)
nonomim_ko$group <- factor(nonomim_ko$group, levels = nonomim_ko$group)
omim_ko_pie <- ggplot(omim_ko, aes(x="", y=value, fill=group))
omim_ko_pie<- omim_ko_pie+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
omim_ko_pie<- omim_ko_pie+scale_fill_manual(values=c("slategray1","slategray"))+theme(legend.position="none")
omim_ko_pie<- omim_ko_pie+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
omim_ko_pie<- omim_ko_pie+ggtitle(paste("OMIM genes \n n=",length(which(universe_df$omim=="Y"))))+theme(plot.title = element_text(size = 20, face = "bold"))

nonomim_ko_pie <- ggplot(nonomim_ko, aes(x="", y=value, fill=group))
nonomim_ko_pie<- nonomim_ko_pie+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
nonomim_ko_pie<- nonomim_ko_pie+scale_fill_manual(values=c("slategray1","slategray"))+theme(legend.position="none")
nonomim_ko_pie<- nonomim_ko_pie+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
nonomim_ko_pie<- nonomim_ko_pie+ggtitle(paste("non-OMIM genes \n n=",length(which(is.na(universe_df$omim)))))+theme(plot.title = element_text(size = 20, face = "bold"))


png("Sandra_Figures/Figs/fig2a-KOpies.png",width=1000,height=500,type="quartz",res=150,bg = "transparent")
ggarrange(omim_ko_pie,nonomim_ko_pie,common.legend=TRUE,legend="bottom")
  dev.off()
  
