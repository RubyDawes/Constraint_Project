load("output/Data/human_lethal_genes.rda")

#Facts
length(which(lethal_genes$mouse_ko=="Y"))
length(which(lethal_genes$lethal_mouse=="Y"))
length(which(lethal_genes$lethal_mouse=="N"))

length(which(universe_df$Inheritance_pattern=="AR"))/length(which(!is.na(universe_df$Inheritance_pattern)))
length(which(lethal_genes$Inheritance_pattern=="AR"))/length(which(!is.na(lethal_genes$Inheritance_pattern)))

length(which(universe_df$lethal_mouse=="Y"))/length(which(!is.na(universe_df$lethal_mouse)))
length(which(lethal_genes$lethal_mouse=="Y"))/length(which(!is.na(lethal_genes$lethal_mouse)))

length(which(lethal_genes$pLI>=0.9))/length(which(!is.na(lethal_genes$pLI)))
length(which(universe_df$pLI>=0.9))/length(which(!is.na(universe_df$pLI)))

length(which(lethal_genes$cell_essential=="Y"))/length(which(!is.na(lethal_genes$cell_essential)))
length(which(universe_df$cell_essential=="Y"))/length(which(!is.na(universe_df$cell_essential)))

#1B: inheritance stacked bar charts
inh_mis <- data.frame(AR=c(length(which(lethal_genes$Inheritance_pattern=="AR"&lethal_genes$mis_z>=3.09))/length(which(lethal_genes$Inheritance_pattern=="AR"&!is.na(lethal_genes$mis_z))),
                           length(which(lethal_genes$Inheritance_pattern=="AR"&lethal_genes$mis_z<3.09))/length(which(lethal_genes$Inheritance_pattern=="AR"&!is.na(lethal_genes$mis_z)))),
                      "ARAD"=c(length(which(lethal_genes$Inheritance_pattern=="AR,AD"&lethal_genes$mis_z>=3.09))/length(which(lethal_genes$Inheritance_pattern=="AR,AD"&!is.na(lethal_genes$mis_z))),
                               length(which(lethal_genes$Inheritance_pattern=="AR,AD"&lethal_genes$mis_z<3.09))/length(which(lethal_genes$Inheritance_pattern=="AR,AD"&!is.na(lethal_genes$mis_z)))),
                      AD=c(length(which(lethal_genes$Inheritance_pattern=="AD"&lethal_genes$mis_z>=3.09))/length(which(lethal_genes$Inheritance_pattern=="AD"&!is.na(lethal_genes$mis_z))),
                           length(which(lethal_genes$Inheritance_pattern=="AD"&lethal_genes$mis_z<3.09))/length(which(lethal_genes$Inheritance_pattern=="AD"&!is.na(lethal_genes$mis_z)))),
                      XL=c(length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&lethal_genes$mis_z>=3.09))/length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&!is.na(lethal_genes$mis_z))),
                           length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&lethal_genes$mis_z<3.09))/length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&!is.na(lethal_genes$mis_z)))))

inh_mism <- melt(inh_mis)
inh_mism$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 4)

f <- ggplot(dat=inh_mism, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="bottom")
ggsave("Sandra_Figures/Figs/fig1bi.pdf",height=7, width=7, units='cm')
rm(inh_mis,inh_mism,f)
inh_lof <- data.frame(AR=c(length(which(lethal_genes$Inheritance_pattern=="AR"&lethal_genes$pLI>=0.9))/length(which(lethal_genes$Inheritance_pattern=="AR"&!is.na(lethal_genes$pLI))),
                           length(which(lethal_genes$Inheritance_pattern=="AR"&lethal_genes$pLI<0.9))/length(which(lethal_genes$Inheritance_pattern=="AR"&!is.na(lethal_genes$pLI)))),
                      "ARAD"=c(length(which(lethal_genes$Inheritance_pattern=="AR,AD"&lethal_genes$pLI>=0.9))/length(which(lethal_genes$Inheritance_pattern=="AR,AD"&!is.na(lethal_genes$pLI))),
                               length(which(lethal_genes$Inheritance_pattern=="AR,AD"&lethal_genes$pLI<0.9))/length(which(lethal_genes$Inheritance_pattern=="AR,AD"&!is.na(lethal_genes$pLI)))),
                      AD=c(length(which(lethal_genes$Inheritance_pattern=="AD"&lethal_genes$pLI>=0.9))/length(which(lethal_genes$Inheritance_pattern=="AD"&!is.na(lethal_genes$pLI))),
                           length(which(lethal_genes$Inheritance_pattern=="AD"&lethal_genes$pLI<0.9))/length(which(lethal_genes$Inheritance_pattern=="AD"&!is.na(lethal_genes$pLI)))),
                      XL=c(length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&lethal_genes$pLI>=0.9))/length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&!is.na(lethal_genes$pLI))),
                           length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&lethal_genes$pLI<0.9))/length(which(grepl("XL",lethal_genes$Inheritance_pattern)&!grepl("MT",lethal_genes$Inheritance_pattern)&!is.na(lethal_genes$pLI)))))


inh_lofm <- melt(inh_lof)
inh_lofm$lof_constraint <- rep(c("LoF constraint     ", "no LoF constraint     "), 4)

e <- ggplot(dat=inh_lofm, aes(x=variable, y=value, fill=lof_constraint))
e<- e+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
e<- e+labs(y = "Proportion")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="bottom")


ggsave("Sandra_Figures/Figs/fig1bii.pdf",height=7, width=7, units='cm')
rm(inh_lof,inh_lofm,e)

#scatter plot
lethal_scatter <- ggplot(lethal_genes[which(!is.na(lethal_genes$mis_z)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,alpha=0.3)+scatter_theme_goodsizes()+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("Non-OMIM genes \n n=",length(which(is.na(universe_df$omim)&!is.na(universe_df$exac))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01))+
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14))
ggsave("Sandra_Figures/Figs/fig1ai.pdf",height=8.9, width=8.9, units='cm')

#venn human vs mouse lethal
universe_df$human_lethal <- ifelse(is.na(vlookup(universe_df$gene,lethal_genes)),"N","Y")

mouse=textGrob("Mouse Lethal", gp=gpar(fontsize=40,fontfamily="Helvetica"),rot=90)
human=textGrob("Human Lethal", gp=gpar(fontsize=40,fontfamily="Helvetica"),rot=270)
a3=draw.pairwise.venn(area1 = length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                      area2 = length(which(universe_df$human_lethal=="Y"&!is.na(universe_df$mouse_ko))), 
                      cross.area=length(which(universe_df$human_lethal=="Y"&universe_df$lethal_mouse=="Y")), 
                      category = c("", ""), lty = "blank",fill = c("midnightblue", "orangered4"), euler.d = TRUE, 
                      scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica")

png('Sandra_Figures/Figs/mouse_vs_human_lethal.png',width=30,height=30,units="cm",res=1000)
grid.arrange(gTree(children=a3),left=mouse,right=human)
dev.off()


