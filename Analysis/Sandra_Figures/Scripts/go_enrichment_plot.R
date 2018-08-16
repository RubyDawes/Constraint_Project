library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)

#all genes go slim annotations- ALL
all <- unlist(universe_df$go_terms[which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$cell_essential))])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,all_count = aa$Count,all_percent=aa$Percent)
rm(aa,all)

#OMIM genes that are either mouse+cell essential, or just mouse lethal (for comparison)
myIds <- unlist(universe_df$go_terms[which(universe_df$omim=="Y"&((universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")|(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")))])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshady$comp_count <- a$Count
slimshady$comp_percent <- a$Percent
rm(a)


#mouse lethal cell unessential OMIM go annotations- CATA
myIds <- unlist(universe_df$go_terms[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")])
myCollection <- GOCollection(myIds)
b<-goSlim(myCollection, slim, "BP")

slimshady$cata_count <- b$Count
slimshady$cata_percent <- b$Percent
rm(b)
#mouse lethal+cell essential OMIM go annotations- CATB
myIds <- unlist(universe_df$go_terms[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")])
myCollection <- GOCollection(myIds)
c<-goSlim(myCollection, slim, "BP")

slimshady$catb_count <- c$Count
slimshady$catb_percent <- c$Percent
rm(c)

#mouse+cell unessential OMIM go annotations- CATC
myIds <- unlist(universe_df$go_terms[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&universe_df$cell_essential=="N")])
myCollection <- GOCollection(myIds)
d<-goSlim(myCollection, slim, "BP")

slimshady$catc_count <- d$Count
slimshady$catc_percent <- d$Percent
rm(d)

slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$cata_percent[which(slimshady$cata_count==0)]<-slimshady$cata_percent[which(slimshady$cata_count==1)][1]
slimshady$catb_percent[which(slimshady$catb_percent==0)]<-slimshady$catb_percent[which(slimshady$catb_count==1)][1]
slimshady$cata_count[which(slimshady$cata_count==0)]<-1
slimshady$catb_count[which(slimshady$catb_count==0)]<-1

slimshady$catacatb_diff <- slimshady$cata_percent-slimshady$catb_percent
slimshady <- slimshady[order(slimshady$catacatb_diff),]


shady_heatmap <- data.frame(go_terms = slimshady$go_term)
shady_heatmap$odds_ratio <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$cata_count[x],(slimshady$cata_count[x]/(slimshady$cata_percent[x]/100))*(1-slimshady$cata_percent[x]/100),
                             slimshady$catb_count[x],(slimshady$catb_count[x]/(slimshady$catb_percent[x]/100))*(1-slimshady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
shady_heatmap$confint_lower <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$cata_count[x],(slimshady$cata_count[x]/(slimshady$cata_percent[x]/100))*(1-slimshady$cata_percent[x]/100),
                             slimshady$catb_count[x],(slimshady$catb_count[x]/(slimshady$catb_percent[x]/100))*(1-slimshady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
shady_heatmap$confint_higher <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$cata_count[x],(slimshady$cata_count[x]/(slimshady$cata_percent[x]/100))*(1-slimshady$cata_percent[x]/100),
                             slimshady$catb_count[x],(slimshady$catb_count[x]/(slimshady$catb_percent[x]/100))*(1-slimshady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
shady_heatmap$pval <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$cata_count[x],(slimshady$cata_count[x]/(slimshady$cata_percent[x]/100))*(1-slimshady$cata_percent[x]/100),
                             slimshady$catb_count[x],(slimshady$catb_count[x]/(slimshady$catb_percent[x]/100))*(1-slimshady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

shady_heatmap <- shady_heatmap[order(shady_heatmap$odds_ratio),]

shady_heatmap$gene <- rep("gene",length(shady_heatmap$go_terms))
shady_heatmap$logodds <- log(shady_heatmap$odds_ratio)
shady_heatmap <- shady_heatmap[-which(shady_heatmap$go_terms=="biological_process"),]
shady_heatmap <- shady_heatmap[-which(shady_heatmap$pval>0.05),]

rm(myCollection,slim,fl,myIds)

#fixing abbreviated GO names
shady_heatmap$go_terms<-as.character(shady_heatmap$go_terms)
shady_heatmap$go_terms[5]<-"nucleobase-containing compound catabolic process"
shady_heatmap$go_terms[15]<-"generation of precursor metabolites and energy"
shady_heatmap$go_terms[18]<-"cellular nitrogen compound metabolic process"
shady_heatmap$go_terms[19]<-"cellular amino acid metabolic process"
shady_heatmap$go_terms[33]<-"anatomical structure formation involved in morphogenesis"
shady_heatmap$go_terms <- factor(shady_heatmap$go_terms, levels = shady_heatmap$go_terms)
#plotting heat map
ggplot(data = shady_heatmap, aes(x= gene,y = go_terms)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "sandybrown",mid="white",high = "firebrick4",trans="log",breaks=c(min(shady_heatmap$odds_ratio),1,max(shady_heatmap$odds_ratio)),labels=c("cell essential","intermediate","mouse lethal"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/go_heatmap.pdf",height=18, width=18, units='cm')


#doing bar graph instead
ggplot(shady_heatmap) +
  geom_bar( aes(x=go_terms, y=odds_ratio,fill=odds_ratio ), stat="identity") +
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10")+scale_fill_gradient2(low = "sandybrown", mid = "white",
                                                         high = "firebrick4",trans="log")+  
  theme(axis.text.x = element_text(angle=60, hjust=1))+coord_flip()+
  ylab("Odds Ratio \n (Cell essential OMIM vs Mouse lethal OMIM)")
ggsave("Analysis/Sandra_Figures/Figs/go_bar.pdf",height=18, width=18, units='cm')
