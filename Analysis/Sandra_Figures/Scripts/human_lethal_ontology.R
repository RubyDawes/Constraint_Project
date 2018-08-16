library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
#slim <- getOBOCollection(fl)

#all oMIM genes go slim annotations- COMB
comb <- unlist(universe_df$go_terms[which(universe_df$human_lethal=="Y"|is.na(universe_df$omim))])
myCollection <- GOCollection(comb)
a<-goSlim(myCollection, slim, "BP")

lethal_shady <- data.frame(go_id=rownames(a),go_term=a$Term,
                        comb_count = a$Count,
                        comb_percent=a$Percent)
rm(a,comb)

#human lethal oMIM genes go slim annotations- CATA
cata <- unlist(universe_df$go_terms[which(universe_df$omim=="Y"&universe_df$human_lethal=="Y")])
myCollection <- GOCollection(cata)
a<-goSlim(myCollection, slim, "BP")
lethal_shady$cata_count <- a$Count
lethal_shady$cata_percent <- a$Percent

# non human lethal oMIM genes go slim annotations- CATB
catb <- unlist(universe_df$go_terms[which(is.na(universe_df$omim))])
myCollection <- GOCollection(catb)
b<-goSlim(myCollection, slim, "BP")
lethal_shady$catb_count <- b$Count
lethal_shady$catb_percent <- b$Percent
rm(a,b,cata,catb)

lethal_shady<-lethal_shady[-c(38,47,69),]

#making 0 counts 1 to avoid infinity OR
lethal_shady$cata_percent[which(lethal_shady$cata_count==0)]<-lethal_shady$cata_percent[which(lethal_shady$cata_count==3)][1]/3
lethal_shady$catb_percent[which(lethal_shady$catb_percent==0)]<-lethal_shady$catb_percent[which(lethal_shady$catb_count==23)][1]/23
lethal_shady$cata_count[which(lethal_shady$cata_count==0)]<-1
lethal_shady$catb_count[which(lethal_shady$catb_count==0)]<-1

lethal_shady$catacatb_diff <- lethal_shady$cata_percent-lethal_shady$catb_percent

lethal_shady$odds_ratio <- unlist(lapply(1:length(lethal_shady$go_term),function(x) {
  fish<-fisher.test(matrix(c(lethal_shady$cata_count[x],(lethal_shady$cata_count[x]/(lethal_shady$cata_percent[x]/100))*(1-lethal_shady$cata_percent[x]/100),
                             lethal_shady$catb_count[x],(lethal_shady$catb_count[x]/(lethal_shady$catb_percent[x]/100))*(1-lethal_shady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))

lethal_shady$confint_lower <- unlist(lapply(1:length(lethal_shady$go_term),function(x) {
  fish<-fisher.test(matrix(c(lethal_shady$cata_count[x],(lethal_shady$cata_count[x]/(lethal_shady$cata_percent[x]/100))*(1-lethal_shady$cata_percent[x]/100),
                             lethal_shady$catb_count[x],(lethal_shady$catb_count[x]/(lethal_shady$catb_percent[x]/100))*(1-lethal_shady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
lethal_shady$confint_higher <- unlist(lapply(1:length(lethal_shady$go_term),function(x) {
  fish<-fisher.test(matrix(c(lethal_shady$cata_count[x],(lethal_shady$cata_count[x]/(lethal_shady$cata_percent[x]/100))*(1-lethal_shady$cata_percent[x]/100),
                             lethal_shady$catb_count[x],(lethal_shady$catb_count[x]/(lethal_shady$catb_percent[x]/100))*(1-lethal_shady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
lethal_shady$pval <- unlist(lapply(1:length(lethal_shady$go_term),function(x) {
  fish<-fisher.test(matrix(c(lethal_shady$cata_count[x],(lethal_shady$cata_count[x]/(lethal_shady$cata_percent[x]/100))*(1-lethal_shady$cata_percent[x]/100),
                             lethal_shady$catb_count[x],(lethal_shady$catb_count[x]/(lethal_shady$catb_percent[x]/100))*(1-lethal_shady$catb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

lethal_shady <- lethal_shady[order(lethal_shady$odds_ratio),]

lethal_shady$gene <- rep("gene",length(lethal_shady$go_term))
lethal_shady$logodds <- log(lethal_shady$odds_ratio)
lethal_shady <- lethal_shady[-which(lethal_shady$go_term=="biological_process"),]
lethal_shady <- lethal_shady[-which(lethal_shady$pval>0.05),]

#fixing abbreviated GO names
lethal_shady$go_term<-as.character(lethal_shady$go_term)
lethal_shady$go_term[8]<-"cellular nitrogen compound metabolic process"
lethal_shady$go_term[12]<-"cellular amino acid metabolic process"
lethal_shady$go_term[18]<-"generation of precursor metabolites and energy"
lethal_shady$go_term <- factor(lethal_shady$go_term, levels = lethal_shady$go_term)

#plotting heat map
ggplot(data = lethal_shady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "sandybrown",mid="white",high = "firebrick4",trans="log",breaks=c(min(lethal_shady$odds_ratio),1,max(lethal_shady$odds_ratio)),labels=c("non-OMIM","intermediate","lethal OMIM"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/go_heatmap.pdf",height=18, width=18, units='cm')

