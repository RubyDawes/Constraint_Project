####loading in curated lethal genes list#####
load("output/Data/human_lethal_genes.rda")
#lethal_mc<-read.xlsx("Gene_lists/lethal_genes/human_lethal_genes_MC.xlsx")

#lethal_genes$listA <- lethal_mc$listA
#lethal_genes$listB <- lethal_mc$listB
#lethal_genes$earliest.death <- lethal_mc$earliest.death
#lethal_genes$caveats <- lethal_mc$caveats
#save(lethal_genes, file="output/Data/human_lethal_genes.rda", compress="bzip2")



####plotting proportion of mouse lethal genes in each list####
mouse_lethal_prop <- data.frame(non_OMIM=c(length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                           length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(is.na(universe_df$omim)&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                           length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(is.na(universe_df$omim)&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen)))),
                                OMIM=c(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                       length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen)))),
                                Human_Lethal_B=c(length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                               length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                               length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_B=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen)))),
                                Human_Lethal_A=c(length(which(universe_df$human_lethal_A=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_A=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                                 length(which(universe_df$human_lethal_A=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_A=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen))),
                                                 length(which(universe_df$human_lethal_A=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_A=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$all_MP_phen)))))

mouse_lethal_propm <- melt(mouse_lethal_prop)
mouse_lethal_propm$mis_constraint <- rep(c("1. Mouse Lethal     ","2. Premature death     ", "3. Mouse non-Lethal     "), 4)

f <- ggplot(dat=mouse_lethal_propm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","grey","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/lethal_prop_MC.pdf",height=15, width=23, units='cm')
rm(mouse_lethal_prop,mouse_lethal_propm,f)



####plotting proportion of constrained genes in each list####
inh_mis_human_lethal_A <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$mis_z>=3.09))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$mis_z<3.09))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   AR=c(length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   "ARAD"=c(length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                            length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   AD=c(length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   XL=c(length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$mis_z>=3.09))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$mis_z<3.09))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))))
inh_mis_human_lethal_Am <- melt(inh_mis_human_lethal_A)
inh_mis_human_lethal_Am$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_mis_human_lethal_Am, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/lethalA_mis_constraint_bar.pdf",height=7, width=7, units='cm')

inh_mis_human_lethal_B <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AR=c(length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     "ARAD"=c(length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                              length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AD=c(length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     XL=c(length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$mis_z>=3.09))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$mis_z<3.09))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))))

inh_mis_human_lethal_Bm <- melt(inh_mis_human_lethal_B)
inh_mis_human_lethal_Bm$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_mis_human_lethal_Bm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="bottom")

ggsave("Analysis/Sandra_Figures/Figs/lethalB_mis_constraint_bar.pdf",height=7, width=7, units='cm')
rm(inh_mis_human_lethal_A,inh_mis_human_lethal_B,inh_mis_human_lethal_Am,inh_mis_human_lethal_Bm,f)

inh_lof_human_lethal_A <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$pLI>=0.9))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$pLI<0.9))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   AR=c(length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   "ARAD"=c(length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                            length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   AD=c(length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                   XL=c(length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$pLI>=0.9))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                        length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&universe_df$pLI<0.9))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))))

inh_lof_human_lethal_Am <- melt(inh_lof_human_lethal_A)
inh_lof_human_lethal_Am$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_lof_human_lethal_Am, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/lethalA_lof_constraint_bar.pdf",height=7, width=7, units='cm')

inh_lof_human_lethal_B <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(grepl(pattern="MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AR=c(length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AR"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     "ARAD"=c(length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                              length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AR,AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     AD=c(length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(universe_df$lethal_inheritance=="AD"&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))),
                                     XL=c(length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI>=0.9))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac))),
                                          length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&universe_df$pLI<0.9))/length(which(grepl("XL",universe_df$lethal_inheritance)&!grepl("MT",universe_df$lethal_inheritance)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$exac)))))

inh_lof_human_lethal_Bm <- melt(inh_lof_human_lethal_B)
inh_lof_human_lethal_Bm$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_lof_human_lethal_Bm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/lethalB_lof_constraint_bar.pdf",height=7, width=7, units='cm')


regconst_lethalA <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_A=="Y"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_B=="Y"&!is.na(universe_df$any_constraint))),
                            length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_A=="Y"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint)))),
                       AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_A=="Y"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint))),
                            length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_A=="Y"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint)))),
                       "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_A=="Y"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint))),
                                length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_A=="Y"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint)))),
                       AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_A=="Y"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint))),
                            length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_A=="Y"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint)))),
                       XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_A=="Y"&(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint))),
                            length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_A=="Y"&!(universe_df$mis_z>=3.09|universe_df$any_reg_constraint=="Y")))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$human_lethal_A=="Y"&!is.na(universe_df$any_constraint)))))


regconst_lethalAm <- melt(regconst_lethalA)
regconst_lethalAm$regconstraint <- rep(c("Regional or whole-gene constraint     ", "no regional or whole-gene constraint     "), 5)

e <- ggplot(dat=regconst_lethalAm, aes(x=variable, y=value, fill=regconstraint))
e<- e+geom_bar(width = 0.8, stat = "identity",color="slategray")+scale_y_continuous(expand = c(0, 0)) +bar_theme()
e<- e+labs(y = "Proportion")+scale_fill_manual(values=c("steelblue3","lightsalmon2"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/lethalA_reg_constraint.pdf",height=7, width=7, units='cm')
rm(regconst,regconstm,e)

####grouping into dominant/recessive####

universe_df$Inheritance_grouped <- ifelse((universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR"),"recessive",
                                          ifelse((universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd"),"dominant","recessive/dominant"))

lethalB_constraint <- data.frame(Recessive=c(length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_grouped=="recessive"&!is.na(universe_df$exac))),
                                            length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_grouped=="recessive"&!is.na(universe_df$exac)))),
                                Mixed=c(length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_grouped=="recessive/dominant"&!is.na(universe_df$exac))),
                                        length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_grouped=="recessive/dominant"&!is.na(universe_df$exac)))),
                                Dominant=c(length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_grouped=="dominant"&!is.na(universe_df$exac))),
                                           length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_grouped=="dominant"&!is.na(universe_df$exac)))))

lethalB_constraintm <- melt(lethalB_constraint)
lethalB_constraintm$constraint <- rep(c("Constraint     ", "No constraint     "), 3)

f <- ggplot(dat=lethalB_constraintm, aes(x=variable, y=value, fill=constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("#CD5555","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/grouped_constraint_lethalB.pdf",height=7, width=7, units='cm')

lethalA_constraint <- data.frame(Recessive=c(length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$human_lethal_A=="Y"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                             length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$human_lethal_A=="Y"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                 Mixed=c(length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$human_lethal_A=="Y"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                         length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$human_lethal_A=="Y"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))),
                                 Dominant=c(length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$human_lethal_A=="Y"&universe_df$constrained=="Y"))/length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac))),
                                            length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$human_lethal_A=="Y"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$human_lethal_A=="Y"&!is.na(universe_df$exac)))))

lethalA_constraintm <- melt(lethalA_constraint)
lethalA_constraintm$constraint <- rep(c("Constraint     ", "No constraint     "), 3)

f <- ggplot(dat=lethalA_constraintm, aes(x=variable, y=value, fill=constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("#CD5555","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/grouped_constraint_lethalA.pdf",height=7, width=7, units='cm')


lethalB_constraint_reg <- data.frame(Recessive=c(length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$any_constraint=="Y"))/length(which(universe_df$Inheritance_grouped=="recessive"&!is.na(universe_df$exac))),
                                             length(which(universe_df$Inheritance_grouped=="recessive"&universe_df$any_constraint=="N"))/length(which(universe_df$Inheritance_grouped=="recessive"&!is.na(universe_df$exac)))),
                                 Mixed=c(length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$any_constraint=="Y"))/length(which(universe_df$Inheritance_grouped=="recessive/dominant"&!is.na(universe_df$exac))),
                                         length(which(universe_df$Inheritance_grouped=="recessive/dominant"&universe_df$any_constraint=="N"))/length(which(universe_df$Inheritance_grouped=="recessive/dominant"&!is.na(universe_df$exac)))),
                                 Dominant=c(length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$any_constraint=="Y"))/length(which(universe_df$Inheritance_grouped=="dominant"&!is.na(universe_df$exac))),
                                            length(which(universe_df$Inheritance_grouped=="dominant"&universe_df$any_constraint=="N"))/length(which(universe_df$Inheritance_grouped=="dominant"&!is.na(universe_df$exac)))))

lethalB_constraint_regm <- melt(lethalB_constraint_reg)
lethalB_constraint_regm$constraint <- rep(c("Constraint     ", "No constraint     "), 3)

f <- ggplot(dat=lethalB_constraint_regm, aes(x=variable, y=value, fill=constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("#CD5555","steelblue3"))+theme(legend.position="bottom")
ggsave("Analysis/Sandra_Figures/Figs/grouped_constraintreg_lethalB.pdf",height=7, width=7, units='cm')




###heatmaps of dominant vs recessive lethal genes GO terms####
library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############lethal genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_A=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshadyBP <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_A=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshadyBP$rec_count <- a$Count
slimshadyBP$rec_percent <- a$Percent
rm(a)

slimshadyBP<-slimshadyBP[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyBP$dom_percent[which(slimshadyBP$dom_count==0)]<-slimshadyBP$dom_percent[which(slimshadyBP$dom_count==1)][1]
slimshadyBP$rec_percent[which(slimshadyBP$rec_percent==0)]<-slimshadyBP$rec_percent[which(slimshadyBP$rec_count==1)][1]
slimshadyBP$dom_count[which(slimshadyBP$dom_count==0)]<-1
slimshadyBP$rec_count[which(slimshadyBP$rec_count==0)]<-1

slimshadyBP$domrec_diff <- slimshadyBP$dom_percent-slimshadyBP$rec_percent


slimshadyBP$odds_ratio <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyBP$confint_lower <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyBP$confint_higher <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyBP$pval <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyBP <- slimshadyBP[-which(slimshadyBP$go_term=="biological_process"),]

slimshadyBP <- slimshadyBP[order(slimshadyBP$odds_ratio),]

slimshadyBP$gene <- rep("gene",length(slimshadyBP$go_term))
#fixing abbreviated GO names
slimshadyBP$go_term<-as.character(slimshadyBP$go_term)
slimshadyBP$go_term[which(slimshadyBP$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshadyBP$go_term[which(slimshadyBP$go_id=="GO:0034655")]<-"nucleobase-containing compound catabolic process"
#slimshadyBP$go_term[which(slimshadyBP$go_id=="GO:0030705")]<-"cytoskeleton-dependent intracellular transport"
#slimshadyBP$go_term[which(slimshadyBP$go_id=="GO:0006091")]<-"generation of precursor metabolites and energy"


slimshadyBP$go_term <- factor(slimshadyBP$go_term, levels = slimshadyBP$go_term)

# SETTING SIGNIFICANCE CUTOFF
slimshadyBP05<-slimshadyBP[which(slimshadyBP$pval<0.05),]
slimshadyBP005<-slimshadyBP[which(slimshadyBP$pval<0.005),]
slimshadyBP0005<-slimshadyBP[which(slimshadyBP$pval<0.0005),]

#plotting heat map
ggplot(data = slimshadyBP05, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyBP$odds_ratio),1,max(slimshadyBP$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_BP_MC_0.05.pdf",height=18, width=18, units='cm')
ggplot(data = slimshadyBP005, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyBP$odds_ratio),1,max(slimshadyBP$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_BP_MC_0.005.pdf",height=18, width=18, units='cm')
ggplot(data = slimshadyBP0005, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyBP$odds_ratio),1,max(slimshadyBP$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_BP_MC_0.0005.pdf",height=18, width=18, units='cm')






############lethal genes dominant vs recessive MF-listB################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "MF")

slimshadyMF <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_B=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "MF")

slimshadyMF$rec_count <- a$Count
slimshadyMF$rec_percent <- a$Percent
rm(a)

slimshadyMF<-slimshadyMF[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyMF$dom_percent[which(slimshadyMF$dom_count==0)]<-slimshadyMF$dom_percent[which(slimshadyMF$dom_count==1)][1]
slimshadyMF$rec_percent[which(slimshadyMF$rec_percent==0)]<-slimshadyMF$rec_percent[which(slimshadyMF$rec_count==1)][1]
slimshadyMF$dom_count[which(slimshadyMF$dom_count==0)]<-1
slimshadyMF$rec_count[which(slimshadyMF$rec_count==0)]<-1

slimshadyMF$domrec_diff <- slimshadyMF$dom_percent-slimshadyMF$rec_percent


slimshadyMF$odds_ratio <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyMF$confint_lower <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyMF$confint_higher <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyMF$pval <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyMF <- slimshadyMF[-which(slimshadyMF$go_term=="molecular_function"),]

slimshadyMF <- slimshadyMF[order(slimshadyMF$odds_ratio),]

slimshadyMF$gene <- rep("gene",length(slimshadyMF$go_term))
#fixing abbreviated GO names
slimshadyMF$go_term<-as.character(slimshadyMF$go_term)
slimshadyMF$go_term[which(slimshadyMF$go_term=="DNA binding transcription factor ac...")]<-"DNA binding transcription factor"
slimshadyMF$go_term[which(slimshadyMF$go_id=="GO:0016757")]<-"transferase activity, transferring glycosyl groups"
slimshadyMF$go_term[which(slimshadyMF$go_id=="GO:0016765")]<-"transferase activity, transferring alkyl or aryl (other than methyl) groups"



slimshadyMF$go_term <- factor(slimshadyMF$go_term, levels = slimshadyMF$go_term)


# SETTING SIGNIFICANCE CUTOFF
slimshadyMF05<-slimshadyMF[which(slimshadyMF$pval<0.05),]
slimshadyMF005<-slimshadyMF[which(slimshadyMF$pval<0.005),]
slimshadyMF0005<-slimshadyMF[which(slimshadyMF$pval<0.0005),]

#plotting heat map
ggplot(data = slimshadyMF05, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyMF$odds_ratio),1,max(slimshadyMF$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_MF_MC_0.05.pdf",height=18, width=18, units='cm')
ggplot(data = slimshadyMF005, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyMF$odds_ratio),1,max(slimshadyMF$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_MF_MC_0.005.pdf",height=18, width=18, units='cm')
ggplot(data = slimshadyMF0005, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyMF$odds_ratio),1,max(slimshadyMF$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_MF_MC_0.0005.pdf",height=18, width=18, units='cm')








############lethal genes dominant vs recessive CC################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_A=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "CC")

slimshadyCC <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_A=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "CC")

slimshadyCC$rec_count <- a$Count
slimshadyCC$rec_percent <- a$Percent
rm(a)

slimshadyCC<-slimshadyCC[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyCC$dom_percent[which(slimshadyCC$dom_count==0)]<-slimshadyCC$dom_percent[which(slimshadyCC$dom_count==1)][1]
slimshadyCC$rec_percent[which(slimshadyCC$rec_percent==0)]<-slimshadyCC$rec_percent[which(slimshadyCC$rec_count==1)][1]
slimshadyCC$dom_count[which(slimshadyCC$dom_count==0)]<-1
slimshadyCC$rec_count[which(slimshadyCC$rec_count==0)]<-1

slimshadyCC$domrec_diff <- slimshadyCC$dom_percent-slimshadyCC$rec_percent


slimshadyCC$odds_ratio <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyCC$confint_lower <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyCC$confint_higher <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyCC$pval <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyCC <- slimshadyCC[-which(slimshadyCC$go_term=="cellular_component"),]

slimshadyCC <- slimshadyCC[order(slimshadyCC$odds_ratio),]

slimshadyCC$gene <- rep("gene",length(slimshadyCC$go_term))

slimshadyCC$go_term <- factor(slimshadyCC$go_term, levels = slimshadyCC$go_term)

# SETTING SIGNIFICANCE CUTOFF
slimshadyCC05<-slimshadyCC[which(slimshadyCC$pval<0.05),]
slimshadyCC005<-slimshadyCC[which(slimshadyCC$pval<0.005),]
slimshadyCC0005<-slimshadyCC[which(slimshadyCC$pval<0.0005),]

#plotting heat map
ggplot(data = slimshadyCC05, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyCC$odds_ratio),1,max(slimshadyCC$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_CC_MC_0.05.pdf",height=18, width=18, units='cm')
ggplot(data = slimshadyCC005, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyCC$odds_ratio),1,max(slimshadyCC$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_CC_MC_0.005.pdf",height=18, width=18, units='cm')
ggplot(data = slimshadyCC0005, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshadyCC$odds_ratio),1,max(slimshadyCC$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_CC_MC_0.0005.pdf",height=18, width=18, units='cm')






