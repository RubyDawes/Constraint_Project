high_phens <- c("adipose tissue phenotype","behavior/neurological phenotype","cardiovascular system phenotype","cellular phenotype",
                "craniofacial phenotype","digestive/alimentary phenotype","embryo phenotype","endocrine/exocrine gland phenotype",
                "growth/size/body region phenotype","hearing/vestibular/ear phenotype","hematopoietic system phenotype","homeostasis/metabolism phenotype",
                "immune system phenotype","integument phenotype","limbs/digits/tail phenotype","liver/biliary system phenotype","mortality/aging",
                "muscle phenotype","neoplasm","nervous system phenotype","normal phenotype","pigmentation phenotype","renal/urinary system phenotype",
                "reproductive system phenotype","respiratory system phenotype","skeleton phenotype","taste/olfaction phenotype","vision/eye phenotype"
)

mp_shady <- data.frame(mp_names = high_phens)

mp_cat_count <- function(mp_name,list_of_phens) {
  return(sum(unlist(lapply(list_of_phens, function(x) length(which(mp_name%in%x=="TRUE"))))))
}
#combination for comparison/ checking
comb<-universe_df$high_MP_phen[which(universe_df$omim=="Y"&((universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")|(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")))]
mp_shady$comb_count <- unlist(lapply(mp_shady$mp_names, function(x) mp_cat_count(x,comb)))

#mouse lethal cell unessential OMIM hpo cat annotations- CATA
cata<-universe_df$high_MP_phen[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")]
mp_shady$cata_count <- unlist(lapply(mp_shady$mp_names, function(x) mp_cat_count(x,cata)))
mp_shady$cata_perc <- mp_shady$cata_count/length(cata)

#mouse lethal+cell essential OMIM go annotations- CATB
catb <- universe_df$high_MP_phen[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")]
mp_shady$catb_count <- unlist(lapply(mp_shady$mp_names, function(x) mp_cat_count(x,catb)))
mp_shady$catb_perc <- mp_shady$catb_count/length(catb)

mp_shady$cat_diff_perc <- mp_shady$cata_perc-mp_shady$catb_perc

mp_shady$odds_ratio <- unlist(lapply(1:length(mp_shady$mp_names),function(x) {
  fish<-fisher.test(matrix(c(mp_shady$cata_count[x],length(cata)-mp_shady$cata_count[x],
                             mp_shady$catb_count[x],length(catb)-mp_shady$catb_count[x]),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))

mp_shady$confint_lower <- unlist(lapply(1:length(mp_shady$mp_names),function(x) {
  fish<-fisher.test(matrix(c(mp_shady$cata_count[x],length(cata)-mp_shady$cata_count[x],
                             mp_shady$catb_count[x],length(catb)-mp_shady$catb_count[x]),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
mp_shady$confint_higher <- unlist(lapply(1:length(mp_shady$mp_names),function(x) {
  fish<-fisher.test(matrix(c(mp_shady$cata_count[x],length(cata)-mp_shady$cata_count[x],
                             mp_shady$catb_count[x],length(catb)-mp_shady$catb_count[x]),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
mp_shady$pval <- unlist(lapply(1:length(mp_shady$mp_names),function(x) {
  fish<-fisher.test(matrix(c(mp_shady$cata_count[x],length(cata)-mp_shady$cata_count[x],
                             mp_shady$catb_count[x],length(catb)-mp_shady$catb_count[x]),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

mp_shady <- mp_shady[order(mp_shady$odds_ratio),]

mp_shady$gene <- rep("gene",length(mp_shady$mp_names))
mp_shady$logodds <- log(mp_shady$odds_ratio)
mp_shady <- mp_shady[-which(mp_shady$mp_names=="mortality/aging"),]
mp_shady <- mp_shady[-which(mp_shady$pval>0.05),]

rm(high_phens,hpo_id,comb,cata,catb)
mp_shady$mp_names <- factor(mp_shady$mp_names, levels = mp_shady$mp_names)

#plotting heat map
ggplot(data = mp_shady, aes(x= gene,y = mp_names)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "sandybrown",mid="white",
                       high = "firebrick4",trans="log",
                       breaks=c(0.465,1,max(mp_shady$odds_ratio)),
                       labels=c("cell essential","intermediate","mouse lethal"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/mp_heatmap.pdf",height=18, width=15, units='cm')

#doing bar graph instead
ggplot(mp_shady) +
  geom_bar( aes(x=mp_names, y=odds_ratio,fill=odds_ratio ), stat="identity") +
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10")+scale_fill_gradient2(low = "sandybrown", mid = "white",
                                                         high = "firebrick4",trans="log")+  
  theme(axis.text.x = element_text(angle=60, hjust=1))+coord_flip()+
  ylab("Odds Ratio \n (Cell essential OMIM vs Mouse lethal OMIM)")
ggsave("Analysis/Sandra_Figures/Figs/mp_bar.pdf",height=18, width=18, units='cm')


