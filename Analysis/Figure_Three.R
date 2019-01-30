source("plot_functions.R")
#3A:OMIM prenatal/infantile lethality genes are strongly associated with lethal murine phenotypes####

mouse_lethal_prop <- data.frame(non_OMIM=c(length(which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y")&universe_df$mouse_ko=="Y"))/
                                             length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse))&grepl("premature death",universe_df$all_MP_phen)))/
                                                    length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse))&!grepl("premature death",universe_df$all_MP_phen)))/
                                             length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))),
                                OMIM=c(length(which(universe_df$omim=="Y"&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y")&universe_df$mouse_ko=="Y"))/
                                         length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y")),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse))&grepl("premature death",universe_df$all_MP_phen)))/
                                         length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y")),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse))&!grepl("premature death",universe_df$all_MP_phen)))/
                                         length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))),
                                Human_Lethal_B=c(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y")&universe_df$mouse_ko=="Y"))/
                                                   length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y")),
                                                 length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse))&grepl("premature death",universe_df$all_MP_phen)))/
                                                   length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y")),
                                                 length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse))&!grepl("premature death",universe_df$all_MP_phen)))/
                                                   length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))))


levels<- c(paste0("non-OMIM \n n = ",length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))),
           paste0("OMIM \n n = ",length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))),
           paste0("Human Lethal \n List B \n n = ",length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))))

mouse_lethal_propm <- melt(mouse_lethal_prop)
mouse_lethal_propm$mis_constraint <- rep(c("Mouse Lethal     ","Premature death     ", "Mouse non-Lethal     "), 3)
mouse_lethal_propm$mis_constraint<-factor(mouse_lethal_propm$mis_constraint,levels=c("Mouse Lethal     ","Premature death     ", "Mouse non-Lethal     "))
mouse_lethal_propm$variable <- rep(levels,each=3)
mouse_lethal_propm$variable <- factor(mouse_lethal_propm$variable, levels = levels)

labels <- rep("",9)
labels[c(1,2,4,5,7,8)]<-c(paste0(round(mouse_lethal_propm[1,2]*100,0),"%"),
                          paste0(round(mouse_lethal_propm[2,2]*100,0),"%"),
                                       paste0(round(mouse_lethal_propm[4,2]*100,0),"%"),
                                                    paste0(round(mouse_lethal_propm[5,2]*100,0),"%"),
                                                                 paste0(round(mouse_lethal_propm[7,2]*100,0),"%"),
                                                                              paste0(round(mouse_lethal_propm[8,2]*100,0),"%"))


f <- ggplot(dat=mouse_lethal_propm, aes(x=variable, y=value, fill=mis_constraint))+
  geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse = TRUE))+
  scale_y_continuous(expand = c(0, 0),limits=c(0,1)) +bar_theme()+
  labs(y = "Proportion")+scale_fill_manual(values=c("black","grey","steelblue3"))+
  theme(legend.position="right")+geom_text(aes(label = labels),colour="white",size=3,position = position_stack(reverse=TRUE,vjust=1))
f
ggsave("output/Figures/3A.pdf",height=7, width=10, units='cm')
rm(mouse_lethal_prop,mouse_lethal_propm,f,levels,labels)

#3B: Flow chart linking disease involvement with genetic constraint and murine non-viable phenotypes####
flowchart <- data.frame(OMIM = c("Y","","","","","","","N","","","","","",""),
                        number= c(length(which(universe_df$omim=="Y")),"","","","","","",
                                 length(which(is.na(universe_df$omim))),"","","","","",""),
                        constraint = c("Y","","","N","","","no data","Y","","","N","","","no data"),
                        number = c(length(which(universe_df$omim=="Y"&universe_df$constrained=="Y")),"","",
                                       length(which(universe_df$omim=="Y"&universe_df$constrained=="N")),"","",
                                       length(which(universe_df$omim=="Y"&is.na(universe_df$exac))),
                                       length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),"","",
                                       length(which(is.na(universe_df$omim)&universe_df$constrained=="N")),"","",
                                       length(which(is.na(universe_df$omim)&is.na(universe_df$exac)))),
                        "mouse lethal" = c("Y","N","no data","Y","N","no data","","Y","N","no data",
                                "Y","N","no data",""),
                        number = c(length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y")&universe_df$mouse_ko=="Y")),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse)))),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&is.na(universe_df$lethal_mouse))),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y")&universe_df$mouse_ko=="Y")),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse)))),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&is.na(universe_df$lethal_mouse))),"",
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y")&universe_df$mouse_ko=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse)))),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&is.na(universe_df$lethal_mouse))),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y")&universe_df$mouse_ko=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&universe_df$lethal_mouse=="N"&(universe_df$lethal_het_mouse=="N"|is.na(universe_df$lethal_het_mouse)))),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&is.na(universe_df$lethal_mouse))),""))

g <- tableGrob(flowchart)
ggplot2::ggsave('output/Figures/3B.pdf',g,height=12.5, width=14.5, units='cm')

rm(flowchart,g)

