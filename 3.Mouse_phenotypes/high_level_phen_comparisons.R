
universe_mouse <- universe_df[which(universe_df$mouse_ko=="Y"),]
rm(universe_df)
high_phens <- c("adipose tissue phenotype","behavior/neurological phenotype","cardiovascular system phenotype","cellular phenotype",
                "craniofacial phenotype","digestive/alimentary phenotype","embryo phenotype","endocrine/exocrine gland phenotype",
                "growth/size/body region phenotype","hearing/vestibular/ear phenotype","hematopoietic system phenotype","homeostasis/metabolism phenotype",
                "immune system phenotype","integument phenotype","limbs/digits/tail phenotype","liver/biliary system phenotype","mortality/aging",
                "muscle phenotype","neoplasm","nervous system phenotype","normal phenotype","pigmentation phenotype","renal/urinary system phenotype",
                "reproductive system phenotype","respiratory system phenotype","skeleton phenotype","taste/olfaction phenotype","vision/eye phenotype"
                )

for (i in 1:length(high_phens)) {
  universe_mouse[,28+i] <- ifelse(grepl(high_phens[i],universe_mouse$high_MP_phen),1,0)
}
rm(i)
colnames(universe_mouse)[29:ncol(universe_mouse)] <- high_phens

universe_mouse$no_high_phens <- rowSums(universe_mouse[,29:ncol(universe_mouse)])

universe_mouse$cell_essential <- ifelse(universe_mouse$cell_essential_hits>=3,"Y","N")

#violin plots 
ggplot(universe_mouse,aes(factor(constrained),no_high_phens))+geom_violin(aes(fill=factor(constrained)))
ggplot(universe_mouse,aes(factor(omim),no_high_phens))+geom_violin(aes(fill=factor(omim)))+blank_theme()
ggplot(universe_mouse,aes(factor(Inheritance_pattern),no_high_phens))+geom_violin(aes(fill=factor(Inheritance_pattern)))
ggplot(universe_mouse,aes(factor(nmd),no_high_phens))+geom_violin(aes(fill=factor(nmd)))
ggplot(universe_mouse,aes(factor(lethal_mouse),no_high_phens))+geom_violin(aes(fill=factor(lethal_mouse)))
ggplot(universe_mouse,aes(factor(cell_essential),no_high_phens))+geom_violin(aes(fill=factor(cell_essential)))

universe_mouse$lof_intol <- ifelse(universe_mouse$pLI>=0.9,"lof intolerant",ifelse(universe_mouse$pLI<=0.1,"lof tolerant","inbetween"))
universe_mouse$lof_intol <- ifelse(universe_mouse$pLI>=0.9,"lof intolerant","lof tolerant")
ggplot(universe_mouse,aes(factor(lof_intol),no_high_phens))+geom_violin(aes(fill=factor(lof_intol)))

high_phens <- c("adipose tissue ","behavior/neurological ","cardiovascular system ","cellular ",
                "craniofacial ","digestive/alimentary ","embryo ","endocrine/exocrine gland ",
                "growth/size/body region ","hearing/vestibular/ear ","hematopoietic system ","homeostasis/metabolism ",
                "immune system ","integument ","limbs/digits/tail ","liver/biliary system ","mortality/aging",
                "muscle ","neoplasm","nervous system ","normal ","pigmentation ","renal/urinary system ",
                "reproductive system ","respiratory system ","skeleton ","taste/olfaction ","vision/eye "
)
#for genes with a certain high_phen, how many were cell essential
percs_cell_ess <- c()
percs_cell_uness <- c()
p_value <- c()
for (i in 1:length(high_phens)){
percs_cell_ess[i] <- length(which(universe_mouse$cell_essential=="Y"&universe_mouse[,28+i]==1))/length(which(universe_mouse$cell_essential=="Y"))
percs_cell_uness[i] <- length(which(universe_mouse$cell_essential=="N"&universe_mouse[,28+i]==1))/length(which(universe_mouse$cell_essential=="N"))
p_value[i] <- fisher.test(matrix(c(length(which(universe_mouse$cell_essential=="Y"&universe_mouse[,28+i]==1)),length(which(universe_mouse$cell_essential=="Y"&universe_mouse[,28+i]==0)),
                                   length(which(universe_mouse$cell_essential=="N"&universe_mouse[,28+i]==1)),length(which(universe_mouse$cell_essential=="N"&universe_mouse[,28+i]==0))),nrow=2,ncol=2),alternative="two.sided")$p.value
}
rm(i)
perc_diff <- percs_cell_ess - percs_cell_uness
a<- data.frame(high_phens,percs_cell_ess,percs_cell_uness,perc_diff,p_value)
ggplot(melt(a[which(a$p_value<0.05),c(1,4)]),aes(x=reorder(high_phens,value),value))+geom_bar(aes(fill=variable),position="dodge",stat="identity")+coord_flip()
b <- ggplot(melt(a[,c(1,4)]),aes(x=reorder(high_phens,value),value))+geom_bar(aes(fill=variable),position="dodge",stat="identity")
b<- b+coord_flip()+bar_twosided_theme()+geom_hline(yintercept=0)+labs(y="Enrichment in cell essential genes",x="High-level phenotype")
b <- b+scale_y_continuous(breaks = pretty(a$perc_diff, n = 6))+ggtitle("Mouse phenotypes of Cell essential genes")
png("output/figures/cellessential_phens.png",width=1800,height=1000,type="quartz",res=150,bg = "transparent")
b
dev.off()
am<- a[-which(a$high_phens=="mortality/aging"),]
am <- melt(am[,c(1,4)])
#for genes with a certain high_phen, how many were lethal in a mouse
percs_mouse_leth <- c()
percs_mouse_nonleth <- c()
p_value2 <- c()
for (i in 1:length(high_phens)){
  percs_mouse_leth[i] <- length(which(universe_mouse$lethal_mouse=="Y"&universe_mouse[,28+i]==1))/length(which(universe_mouse$lethal_mouse=="Y"))
  percs_mouse_nonleth[i] <- length(which(universe_mouse$lethal_mouse=="N"&universe_mouse[,28+i]==1))/length(which(universe_mouse$lethal_mouse=="N"))
  p_value2[i] <- fisher.test(matrix(c(length(which(universe_mouse$lethal_mouse=="Y"&universe_mouse[,28+i]==1)),length(which(universe_mouse$lethal_mouse=="Y"&universe_mouse[,28+i]==0)),
                                     length(which(universe_mouse$lethal_mouse=="N"&universe_mouse[,28+i]==1)),length(which(universe_mouse$lethal_mouse=="N"&universe_mouse[,28+i]==0))),nrow=2,ncol=2),alternative="two.sided")$p.value
}
rm(i)
perc_diff2 <- percs_mouse_leth - percs_mouse_nonleth
c<- data.frame(high_phens,percs_mouse_leth,percs_mouse_nonleth,perc_diff2,p_value2)
c<- c[-which(c$high_phens=="mortality/aging"),]
d <- ggplot(melt(c[,c(1,4)]),aes(x=reorder(high_phens,am$value),value))+geom_bar(aes(fill=variable),position="dodge",stat="identity")
d<- d+coord_flip()+bar_twosided_theme()+geom_hline(yintercept=0)+labs(y="Enrichment in mouse lethal genes",x="High-level phenotype")
d <- d+scale_y_continuous(breaks = pretty(a$perc_diff, n = 6))+ggtitle("Mouse phenotypes of Mouse lethal genes")

png("output/figures/mouselethal_phens.png",width=1800,height=1000,type="quartz",res=150,bg = "transparent")
d
dev.off()

#for genes with a certain high_phen, how many were lof intolerant/tolerant
percs_lof_intol <- c()
percs_lof_tol <- c()
p_value3 <- c()
for (i in 1:length(high_phens)){
  percs_lof_intol[i] <- length(which(universe_mouse$pLI>=0.9&universe_mouse[,28+i]==1))/length(which(universe_mouse$pLI>=0.9))
  percs_lof_tol[i] <- length(which(universe_mouse$pLI<=0.1&universe_mouse[,28+i]==1))/length(which(universe_mouse$pLI<=0.1))
  p_value3[i] <- fisher.test(matrix(c(length(which(universe_mouse$pLI>=0.9&universe_mouse[,28+i]==1)),length(which(universe_mouse$pLI>=0.9&universe_mouse[,28+i]==0)),
                                      length(which(universe_mouse$pLI<=0.1&universe_mouse[,28+i]==1)),length(which(universe_mouse$pLI<=0.1&universe_mouse[,28+i]==0))),nrow=2,ncol=2),alternative="two.sided")$p.value
}
rm(i)
perc_diff3 <- percs_lof_intol - percs_lof_tol
e<- data.frame(high_phens,percs_lof_intol,percs_lof_tol,perc_diff3,p_value3)
f <- ggplot(melt(c[,c(1,4)]),aes(x=reorder(high_phens,value),value))+geom_bar(aes(fill=variable),position="dodge",stat="identity")
f<- f+coord_flip()+bar_twosided_theme()+geom_hline(yintercept=0)+labs(y="Enrichment in LoF intolerant genes genes",x="High-level phenotype")
f <- f+scale_y_continuous(breaks = pretty(a$perc_diff, n = 6))+ggtitle("Mouse phenotypes of LoF intolerant genes genes")
png("output/figures/lofintolerant_phens.png",width=1800,height=1000,type="quartz",res=150,bg = "transparent")
f
dev.off()

#cell essential vs mouse lethal
cell_vs_mouse <- percs_cell_ess-percs_mouse_leth
p_valuecellvsmouse <- c()
for (i in 1:length(high_phens)){
p_valuecellvsmouse[i] <- fisher.test(matrix(c(length(which(universe_mouse$cell_essential=="Y"&universe_mouse[,28+i]==1)),length(which(universe_mouse$cell_essential=="Y"&universe_mouse[,28+i]==0)),
                                   length(which(universe_mouse$lethal_mouse=="Y"&universe_mouse[,28+i]==1)),length(which(universe_mouse$lethal_mouse=="N"&universe_mouse[,28+i]==0))),nrow=2,ncol=2),alternative="two.sided")$p.value
}
f <- data.frame(high_phens,percs_cell_ess,percs_mouse_leth,cell_vs_mouse,p_valuecellvsmouse)
f<-f[-which(high_phens=="mortality/aging"),]
g <- ggplot(melt(f[,c(1,4)]),aes(x=reorder(high_phens,value),value))+geom_bar(aes(fill=variable),position="dodge",stat="identity")
g<- g+coord_flip()+bar_twosided_theme()+geom_hline(yintercept=0)+labs(y="Cell vs Mouse proportion",x="High-level phenotype")
g <- g+scale_y_continuous(breaks = pretty(g$perc_diff, n = 6))

#cell essential more likely than mouse to have cancer associated phen? are cell essnetial genes cancer genes?
