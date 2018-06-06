
#how many genes with msisense and lof constraint
#source("2.ExAC_constraint/exac_constraint.R")
source("plot_functions.R")

############Figure 1
constrained <- exac$gene[which(exac$mis_z>=3.09|exac$pLI>=0.9)]
non_constrained <- exac$gene[which(exac$mis_z<3.09&exac$pLI<0.9)]

constrained_momim <- exac$gene[which(exac$disease[which(exac$gene%in%constrained)]=="Y")]
constrained_non_momim <- exac$gene[which(exac$disease[which(exac$gene%in%constrained)]=="N")]

non_constrained_momim <- exac$gene[which(exac$disease[which(exac$gene%in%non_constrained)]=="Y")]
non_constrained_non_momim <- exac$gene[which(exac$disease[which(exac$gene%in%non_constrained)]=="N")]

constraint <- data.frame(group=c("Constrained OMIM","Constrained non-OMIM","Non Constrained non-OMIM","Non Constrained OMIM"),
                         value=c(length(constrained_momim),length(constrained_non_momim),length(non_constrained_non_momim),length(non_constrained_momim)))

constraint$group <- factor(constraint$group, levels = constraint$group)

a <- ggplot(constraint, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_manual(values=alpha(c("brown3","brown1","dodgerblue1","dodgerblue3"),0.7))#+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = value),size=4,position = position_stack(vjust = 0.5))

png("Sandra_Figures/Figs/fig1.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
a
dev.off()

###############Figure 2
omim_AR <- exac[which(exac$gene%in%omim$Gene[which(omim$Inheritance_pattern=="AR")]),]
omim_ARAD <- exac[which(exac$gene%in%omim$Gene[which(omim$Inheritance_pattern=="AR,AD")]),]
omim_AD <- exac[which(exac$gene%in%omim$Gene[which(omim$Inheritance_pattern=="AD")]),]
omim_XL <-exac[which(exac$gene%in%omim$Gene[which(grepl(pattern="XL",omim$Inheritance_pattern,ignore.case=FALSE)&!grepl(pattern="MT",omim$Inheritance_pattern,ignore.case=FALSE))]),]
omim_MT <- exac[which(exac$gene%in%omim$Gene[which(grepl(pattern="MT",omim$Inheritance_pattern,ignore.case=FALSE))]),]

inh_mis <- data.frame(MT=c(length(which(omim_MT$mis_z>=3.09))/length(omim_MT$gene),length(which(omim_MT$mis_z<3.09))/length(omim_MT$gene)),
                      AR=c(length(which(omim_AR$mis_z>=3.09))/length(omim_AR$gene),length(which(omim_AR$mis_z<3.09))/length(omim_AR$gene)),
                      ARAD=c(length(which(omim_ARAD$mis_z>=3.09))/length(omim_ARAD$gene),length(which(omim_ARAD$mis_z<3.09))/length(omim_ARAD$gene)),
                      AD=c(length(which(omim_AD$mis_z>=3.09))/length(omim_AD$gene),length(which(omim_AD$mis_z<3.09))/length(omim_AD$gene)),
                      XL=c(length(which(omim_XL$mis_z>=3.09))/length(omim_XL$gene),length(which(omim_XL$mis_z<3.09))/length(omim_XL$gene)))
inh_mism <- melt(inh_mis)
inh_mism$mis_constraint <- rep(c("missense constraint", "no missense constraint"), 5)

d <- ggplot(dat=inh_mism, aes(x=variable, y=value, fill=mis_constraint))
d<- d+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
d<- d+labs(y = "Proportion")+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="bottom")


png("Sandra_Figures/Figs/fig2a.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
d
dev.off()

inh_lof <- data.frame(MT=c(length(which(omim_MT$pLI>=0.9))/length(omim_MT$gene),length(which(omim_MT$pLI<0.9))/length(omim_MT$gene)),
                      AR=c(length(which(omim_AR$pLI>=0.9))/length(omim_AR$gene),length(which(omim_AR$pLI<0.9))/length(omim_AR$gene)),
                      ARAD=c(length(which(omim_ARAD$pLI>=0.9))/length(omim_ARAD$gene),length(which(omim_ARAD$pLI<0.9))/length(omim_ARAD$gene)),
                      AD=c(length(which(omim_AD$pLI>=0.9))/length(omim_AD$gene),length(which(omim_AD$pLI<0.9))/length(omim_AD$gene)),
                      XL=c(length(which(omim_XL$pLI>=0.9))/length(omim_XL$gene),length(which(omim_XL$pLI<0.9))/length(omim_XL$gene)))
inh_lofm <- melt(inh_lof)
inh_lofm$lof_constraint <- rep(c("LoF constraint", "no LoF constraint"), 5)

e <- ggplot(dat=inh_lofm, aes(x=variable, y=value, fill=lof_constraint))
e<- e+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
e<- e+labs(y = "Proportion")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="bottom")


png("Sandra_Figures/Figs/fig2b.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
e
dev.off()

##############Figure 3
load("universe_df.rda")
boy <- data.frame(group=c("disease genes with no mouse phenotype","disease genes with non-lethal mouse phenotype"
                          ,"disease genes with lethal mouse phenotype","nondisease genes with lethal mouse phenotype"
                          ,"nondisease genes with non-lethal mouse phenotype","nondisease genes with no mouse phenotype"),
                  value=c(length(which(universe_df$omim=="Y"&is.na(universe_df$mouse_ko))),length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"&universe_df$lethal_mouse=="N"))
                          ,length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"&universe_df$lethal_mouse=="Y")),length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"&universe_df$lethal_mouse=="Y"))
                          ,length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"&universe_df$lethal_mouse=="N")),length(which(is.na(universe_df$omim)&is.na(universe_df$mouse_ko)))))

boy$group <- factor(boy$group, levels = boy$group)

b <- ggplot(boy, aes(x="", y=value, fill=group))
b<- b+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
b<- b+scale_fill_manual(values=alpha(c("brown3","brown2","brown1","dodgerblue1","dodgerblue2","dodgerblue3"),0.7))#+theme(legend.position="none")
b<- b+geom_text(aes(y = value,label = value),size=4,position = position_stack(vjust = 0.5))

 

##############Figure 4
universe_df$constrained <- ifelse(universe_df$mis_z>=3.09|universe_df$pLI>=0.9,"Y","N")
#constrained genes that are in omim and lethal in a mouse
length(which(universe_df$constrained=="Y"&universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))
#constrained genes that are in omim and non-lethal in mouse or no phenotype
length(which(universe_df$constrained=="Y"&universe_df$omim=="Y"&universe_df$lethal_mouse=="N"))
length(which(universe_df$constrained=="Y"&universe_df$omim=="Y"&is.na(universe_df$mouse_ko)))
#constrained genes that aren't in omim and lethal in a mouse
length(which(universe_df$constrained=="Y"&is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))
#constrained genes that aren't in omim and non-lethal in mouse or no phenotype
length(which(universe_df$constrained=="Y"&is.na(universe_df$omim)&universe_df$lethal_mouse=="N"))
length(which(universe_df$constrained=="Y"&is.na(universe_df$omim)&is.na(universe_df$mouse_ko)))

papoose <- data.frame(group=c("constrained genes that are in omim and non-lethal in mouse or no phenotype","constrained genes that are in omim and lethal in a mouse"
                              ,"constrained genes that aren't in omim and lethal in a mouse","constrained genes that aren't in omim and non-lethal in mouse or no phenotype"),
                      value=c(length(which(universe_df$constrained=="Y"&universe_df$omim=="Y"&universe_df$lethal_mouse=="N"))+length(which(universe_df$constrained=="Y"&universe_df$omim=="Y"&is.na(universe_df$mouse_ko))),length(which(universe_df$constrained=="Y"&universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))
                              ,length(which(universe_df$constrained=="Y"&is.na(universe_df$omim)&universe_df$lethal_mouse=="N"))+length(which(universe_df$constrained=="Y"&is.na(universe_df$omim)&is.na(universe_df$mouse_ko))),length(which(universe_df$constrained=="Y"&is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))))
papoose$group <- factor(papoose$group, levels = papoose$group)

c <- ggplot(papoose, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=alpha(c("brown3","brown1","dodgerblue1","dodgerblue3"),0.7))#+theme(legend.position="none")
c<- c+geom_text(aes(y = value,label = value),size=4,position = position_stack(vjust = 0.5))

png("Sandra_Figures/Figs/fig4.png",width=1500,height=1000,type="quartz",res=150,bg = "transparent")
c
dev.off()

