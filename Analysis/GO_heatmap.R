load("output/Data/GOSLIM.rda")
#getting odds ratios etc
#omim dom vs rec####
annotated_domomim_no <-length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")),"gene"]%in%goslim_anno_V3$V3))
annotated_recomim_no <-length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))

go_fish_V3$domvsrecomim_odds_ratio <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$omimdom_gene_no[x],go_fish_V3$omimrec_gene_no[x],
                             annotated_domomim_no-go_fish_V3$omimdom_gene_no[x],annotated_recomim_no-go_fish_V3$omimrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
    return(fish$estimate[[1]])
}))

go_fish_V3$domvsrecomim_confint_lower <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$omimdom_gene_no[x],go_fish_V3$omimrec_gene_no[x],
                             annotated_domomim_no-go_fish_V3$omimdom_gene_no[x],annotated_recomim_no-go_fish_V3$omimrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))

go_fish_V3$domvsrecomim_confint_higher <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$omimdom_gene_no[x],go_fish_V3$omimrec_gene_no[x],
                             annotated_domomim_no-go_fish_V3$omimdom_gene_no[x],annotated_recomim_no-go_fish_V3$omimrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))

go_fish_V3$domvsrecomim_pval <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$omimdom_gene_no[x],go_fish_V3$omimrec_gene_no[x],
                             annotated_domomim_no-go_fish_V3$omimdom_gene_no[x],annotated_recomim_no-go_fish_V3$omimrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

#lethal dom vs rec####
annotated_domlethal_no <-length(which(universe_df[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),"gene"]%in%goslim_anno_V3$V3))
annotated_reclethal_no <-length(which(universe_df[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))

go_fish_V3$domvsreclethal_odds_ratio <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$lethaldom_gene_no[x],go_fish_V3$lethalrec_gene_no[x],
                             annotated_domlethal_no-go_fish_V3$lethaldom_gene_no[x],annotated_reclethal_no-go_fish_V3$lethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$estimate[[1]])
}))

go_fish_V3$domvsreclethal_confint_lower <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$lethaldom_gene_no[x],go_fish_V3$lethalrec_gene_no[x],
                             annotated_domlethal_no-go_fish_V3$lethaldom_gene_no[x],annotated_reclethal_no-go_fish_V3$lethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))

go_fish_V3$domvsreclethal_confint_higher <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$lethaldom_gene_no[x],go_fish_V3$lethalrec_gene_no[x],
                             annotated_domlethal_no-go_fish_V3$lethaldom_gene_no[x],annotated_reclethal_no-go_fish_V3$lethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))

go_fish_V3$domvsreclethal_pval <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$lethaldom_gene_no[x],go_fish_V3$lethalrec_gene_no[x],
                             annotated_domlethal_no-go_fish_V3$lethaldom_gene_no[x],annotated_reclethal_no-go_fish_V3$lethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


#put lethal dom vs rec####
annotated_putdomlethal_no <-length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"),"gene"]%in%goslim_anno_V3$V3))
annotated_putreclethal_no <-length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"),"gene"]%in%goslim_anno_V3$V3))

go_fish_V3$putdomvsreclethal_odds_ratio <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$putlethaldom_gene_no[x],go_fish_V3$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish_V3$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish_V3$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$estimate[[1]])
}))

go_fish_V3$putdomvsreclethal_confint_lower <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$putlethaldom_gene_no[x],go_fish_V3$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish_V3$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish_V3$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))

go_fish_V3$putdomvsreclethal_confint_higher <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$putlethaldom_gene_no[x],go_fish_V3$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish_V3$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish_V3$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))

go_fish_V3$putdomvsreclethal_pval <- unlist(lapply(1:length(go_fish_V3$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish_V3$putlethaldom_gene_no[x],go_fish_V3$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish_V3$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish_V3$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

#plotting heatmaps####
go_fish_V3$gene <- rep("gene",length(go_fish_V3$go_id))
min<-min(go_fish_V3$putdomvsreclethal_odds_ratio)
max<-max(go_fish_V3$putdomvsreclethal_odds_ratio)
#omim dom vs rec
go_fish_CC <- go_fish_V3[order(go_fish_V3$domvsrecomim_odds_ratio[1:10]),]
go_fish_CC$go_name <- factor(go_fish_CC$go_name,levels=go_fish_CC$go_name)
a<-ggplot(data = go_fish_CC, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = domvsrecomim_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,3.0,round(max,1)),labels=c(round(min,1),0.6,1.0,3.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))


go_fish_MF <- go_fish_V3[order(go_fish_V3$domvsrecomim_odds_ratio[11:20])+10,]
go_fish_MF$go_name <- factor(go_fish_MF$go_name,levels=go_fish_MF$go_name)
b<-ggplot(data = go_fish_MF, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = domvsrecomim_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,3.0,round(max,1)),labels=c(round(min,1),0.6,1.0,3.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

go_fish_BP <- go_fish_V3[order(go_fish_V3$domvsrecomim_odds_ratio[21:30])+20,]
go_fish_BP$go_name <- factor(go_fish_BP$go_name,levels=go_fish_BP$go_name)
c<-ggplot(data = go_fish_BP, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = domvsrecomim_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,3.0,round(max,1)),labels=c(round(min,1),0.6,1.0,3.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

#lethal dom vs rec
go_fish_lethal_CC <- go_fish_V3[order(go_fish_V3$domvsreclethal_odds_ratio[1:10]),]
go_fish_lethal_CC$go_name <- factor(go_fish_lethal_CC$go_name,levels=go_fish_lethal_CC$go_name)
d<-ggplot(data = go_fish_lethal_CC, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = domvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,3.0,round(max,1)),labels=c(round(min,1),0.6,1.0,3.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position = 'none',axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

go_fish_lethal_MF <- go_fish_V3[order(go_fish_V3$domvsreclethal_odds_ratio[11:20])+10,]
go_fish_lethal_MF$go_name <- factor(go_fish_lethal_MF$go_name,levels=go_fish_lethal_MF$go_name)
e<-ggplot(data = go_fish_lethal_MF, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = domvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,3.0,round(max,1)),labels=c(round(min,1),0.6,1.0,3.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

go_fish_lethal_BP <- go_fish_V3[order(go_fish_V3$domvsreclethal_odds_ratio[21:30])+20,]
go_fish_lethal_BP$go_name <- factor(go_fish_lethal_BP$go_name,levels=go_fish_lethal_BP$go_name)
f<-ggplot(data = go_fish_lethal_BP, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = domvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,3.0,round(max,1)),labels=c(round(min,1),0.6,1.0,3.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=9))+
  scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))


#putative dom vs rec
go_fish_putlethal_CC <- go_fish_V3[order(go_fish_V3$putdomvsreclethal_odds_ratio[1:10]),]
go_fish_putlethal_CC$go_name <- factor(go_fish_putlethal_CC$go_name,levels=go_fish_putlethal_CC$go_name)
g<-ggplot(data = go_fish_putlethal_CC, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = putdomvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,round(max,1)),labels=c(round(min,1),0.6,1.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position = 'none',axis.title.y=element_blank(),axis.text.y=element_text(size=9))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

go_fish_putlethal_MF <- go_fish_V3[order(go_fish_V3$putdomvsreclethal_odds_ratio[11:20])+10,]
go_fish_putlethal_MF$go_name <- factor(go_fish_putlethal_MF$go_name,levels=go_fish_putlethal_MF$go_name)
h<-ggplot(data = go_fish_putlethal_MF[which(go_fish_putlethal_MF$putdomvsreclethal_odds_ratio<0.9|go_fish_putlethal_MF$putdomvsreclethal_odds_ratio>1.1),], aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = putdomvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1),0.6,1.0,round(max,1)),labels=c(round(min,1),0.6,1.0,round(max,1)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position = 'none',axis.title.y=element_blank(),axis.text.y=element_text(size=9))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

go_fish_putlethal_BP <- go_fish_V3[order(go_fish_V3$putdomvsreclethal_odds_ratio[21:30])+20,]
go_fish_putlethal_BP$go_name <- factor(go_fish_putlethal_BP$go_name,levels=go_fish_putlethal_BP$go_name)
i<-ggplot(data = go_fish_putlethal_BP, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = putdomvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradient2(low="#a020f0",mid="white",high="#008b45", 
                       breaks=c(round(min,1)+0.01,0.6,1.0,round(max,0)),labels=c(round(min,1),0.6,1.0,round(max,0)),limits=c(min,max),
                       guide = "colourbar",trans="log")+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=9))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))



ggarrange(g,h,i,widths=c(1,1,1.3),ncol=3)
ggsave("output/Figures/5A.pdf",height=10 , width=16, units='cm')



