
source("2.ExAC_constraint/exac_constraint.R")
source("plot_functions.R")
#how many disease genes have constraint? how many non-disease genes have constraint?

omim_mis_constraint <- data.frame(group=c("mis_z>=3.09","mis_z<3.09"),value=c(length(which(exac$mis_z[which(exac$omim=="Y")]>=3.09)),length(which(exac$mis_z[which(exac$omim=="Y")]<3.09))))
a <- ggplot(omim_mis_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("OMIM genes \n n=",length(which(exac$omim=="Y"))))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
a<- a+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value)),family="Avenir"),size=5,position = position_stack(vjust = 0.5))

nmd_mis_constraint <- data.frame(group=c("mis_z>=3.09","mis_z<3.09"),value=c(length(which(exac$mis_z[which(exac$nmd=="Y")]>=3.09)),length(which(exac$mis_z[which(exac$nmd=="Y")]<3.09))))
b <- ggplot(nmd_mis_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("NMD genes\n n=",length(which(exac$nmd=="Y"))))
b<- b+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
b<- b+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=5,position = position_stack(vjust = 0.5))

disease_mis_constraint <- data.frame(group=c("mis_z>=3.09","mis_z<3.09"),value=c(length(which(exac$mis_z[which(exac$disease=="Y")]>=3.09)),length(which(exac$mis_z[which(exac$disease=="Y")]<3.09))))
c <- ggplot(disease_mis_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("disease genes\n (OMIM + NMD)\n n=",length(which(exac$disease=="Y"))))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
c<- c+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=5,position = position_stack(vjust = 0.5))

non_disease_mis_constraint <- data.frame(group=c("mis_z>=3.09","mis_z<3.09"),value=c(length(which(exac$mis_z[which(exac$disease=="N")]>=3.09)),length(which(exac$mis_z[which(exac$disease=="N")]<3.09))))
d <- ggplot(non_disease_mis_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("non-disease genes\n (not in OMIM or NMD)\n n=",length(which(exac$disease=="N"))))
d<- d+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
d<- d+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))

png("output/figures/fig1.png",width=1000,height=1000,type="quartz",res=150)
grid_arrange_shared_legend(a,b,c,d)
dev.off()



#lof constraint
omim_lof_constraint <- data.frame(group=c("pLI>=0.9","pLI<0.9"),value=c(length(which(exac$pLI[which(exac$omim=="Y")]>=0.9)),length(which(exac$pLI[which(exac$omim=="Y")]<0.9))))
e <- ggplot(omim_lof_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("OMIM genes \n n=",length(which(exac$omim=="Y"))))
e<- e+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
e<- e+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))


nmd_lof_constraint <- data.frame(group=c("pLI>=0.9","pLI<0.9"),value=c(length(which(exac$pLI[which(exac$nmd=="Y")]>=0.9)),length(which(exac$pLI[which(exac$nmd=="Y")]<0.9))))
f <- ggplot(nmd_lof_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("NMD genes\n n=",length(which(exac$nmd=="Y"))))
f<- f+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
f<- f+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))

disease_lof_constraint <- data.frame(group=c("pLI>=0.9","pLI<0.9"),value=c(length(which(exac$pLI[which(exac$disease=="Y")]>=0.9)),length(which(exac$pLI[which(exac$disease=="Y")]<0.9))))
g <- ggplot(disease_lof_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("disease genes\n (OMIM + NMD)\n n=",length(which(exac$disease=="Y"))))
g<- g+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
g<- g+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))


non_disease_lof_constraint <- data.frame(group=c("pLI>=0.9","pLI<0.9"),value=c(length(which(exac$pLI[which(exac$disease=="N")]>=0.9)),length(which(exac$pLI[which(exac$disease=="N")]<0.9))))
h <- ggplot(non_disease_lof_constraint, aes(x="", y=value, fill=group))+ggtitle(paste("non-disease genes\n (not in OMIM or NMD)\n n=",length(which(exac$disease=="N"))))
h<- h+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
h<- h+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))

png("output/figures/fig2.png",width=1000,height=1000,type="quartz",res=150)
grid_arrange_shared_legend(e,f,g,h)
dev.off()

#how many constrained genes are disease genes?
mis_constraint_omim <- data.frame(group=c("OMIM","non-OMIM"),value=c(length(which(exac$omim[which(exac$mis_z>=3.09)]=="Y")),length(which(exac$omim[which(exac$mis_z>=3.09)]=="N"))))
i <- ggplot(mis_constraint_omim, aes(x="", y=value, fill=group))+ggtitle(paste("Genes with missense constraint \n n=",length(which(exac$mis_z>=3.09))))
i<- i+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
i<- i+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))
  
lof_constraint_omim <- data.frame(group=c("OMIM","non-OMIM"),value=c(length(which(exac$omim[which(exac$pLI>=0.9)]=="Y")),length(which(exac$omim[which(exac$pLI>=0.9)]=="N"))))
j <-  ggplot(lof_constraint_omim, aes(x="", y=value, fill=group))+ggtitle(paste("Genes with LoF constraint \n n=",length(which(exac$pLI>=0.9))))
j<- j+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
j<- j+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))

constraint_omim <- data.frame(group=c("OMIM","non-OMIM"),value=c(length(which(exac$omim[which(exac$mis_z>=3.09&exac$pLI>=0.9)]=="Y")),length(which(exac$omim[which(exac$mis_z>=3.09&exac$pLI>=0.9)]=="N"))))
k <-  ggplot(constraint_omim, aes(x="", y=value, fill=group))+ggtitle(paste("Genes with missense \n and LoF constraint \n n=",length(which(exac$mis_z>=3.09&exac$pLI>=0.9))))
k<- k+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
k<- k+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))

no_constraint_omim <- data.frame(group=c("OMIM","non-OMIM"),value=c(length(which(exac$omim[which(exac$mis_z<3.09&exac$pLI<0.9)]=="Y")),length(which(exac$omim[which(exac$mis_z<3.09&exac$pLI<0.9)]=="N"))))
l <- ggplot(no_constraint_omim, aes(x="", y=value, fill=group))+ggtitle(paste("Genes with no \n missense or LoF constraint \n n=",length(which(exac$mis_z<3.09&exac$pLI<0.9))))
l<- l+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar("y", start=0)
l<- l+scale_fill_manual(values=c("dodgerblue3","brown3"))+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))

png("output/figures/fig3.png",width=1000,height=1000,type="quartz",res=150)
grid_arrange_shared_legend(i,j,k,l)
dev.off()

