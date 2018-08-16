library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############lethal genes################
#all genes go slim annotations- human lethal
all <- unlist(universe_df$go_terms[which(universe_df$human_lethal=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,lethal_count = aa$Count,lethal_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations
myIds <- unlist(universe_df$go_terms[which(lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshady$all_count <- a$Count
slimshady$all_percent <- a$Percent
rm(a)


slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$lethal_percent[which(slimshady$lethal_count==0)]<-slimshady$lethal_percent[which(slimshady$lethal_count==3)][1]/3
slimshady$all_percent[which(slimshady$all_percent==0)]<-slimshady$all_percent[which(slimshady$all_count==1)][1]
slimshady$lethal_count[which(slimshady$lethal_count==0)]<-1
slimshady$all_count[which(slimshady$all_count==0)]<-1

slimshady$lethalall_diff <- slimshady$lethal_percent-slimshady$all_percent
slimshady <- slimshady[order(slimshady$lethalall_diff),]


slimshady$odds_ratio <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
slimshady$logodds <- log(slimshady$odds_ratio)
slimshady <- slimshady[-which(slimshady$go_term=="biological_process"),]
slimshady <- slimshady[-which(slimshady$pval>0.05),]

rm(myCollection,slim,fl,myIds)

#fixing abbreviated GO names
slimshady$go_terms<-as.character(slimshady$go_term)
slimshady$go_terms[5]<-"nucleobase-containing compound catabolic process"
slimshady$go_terms[15]<-"generation of precursor metabolites and energy"
slimshady$go_terms[18]<-"cellular nitrogen compound metabolic process"
slimshady$go_terms[19]<-"cellular amino acid metabolic process"
slimshady$go_terms[33]<-"anatomical structure formation involved in morphogenesis"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "sandybrown",mid="white",high = "firebrick4",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("depleted in human lethal genes","intermediate","enriched in human lethal genes"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/go_heatmap.pdf",height=18, width=18, units='cm')


#doing bar graph instead
ggplot(slimshady) +
  geom_bar( aes(x=go_term, y=odds_ratio,fill=odds_ratio ), stat="identity") +
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10")+scale_fill_gradient2(low = "sandybrown", mid = "white",
                                                         high = "firebrick4",trans="log")+  
  theme(axis.text.x = element_text(angle=60, hjust=1))+coord_flip()+
  ylab("Odds Ratio \n (Cell essential OMIM vs Mouse lethal OMIM)")
ggsave("Analysis/Sandra_Figures/Figs/go_bar.pdf",height=18, width=18, units='cm')

############candidate genes################
#all genes go slim annotations- candidate genes
all <- unlist(universe_df$go_terms[which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|universe_df$cell_essential=="Y")
                                         &lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,candidate_count = aa$Count,candidate_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations
myIds <- unlist(universe_df$go_terms[which(lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshady$all_count <- a$Count
slimshady$all_percent <- a$Percent
rm(a)


slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$candidate_percent[which(slimshady$candidate_count==0)]<-slimshady$candidate_percent[which(slimshady$candidate_count==3)][1]/3
slimshady$all_percent[which(slimshady$all_percent==0)]<-slimshady$all_percent[which(slimshady$all_count==1)][1]
slimshady$candidate_count[which(slimshady$candidate_count==0)]<-1
slimshady$all_count[which(slimshady$all_count==0)]<-1

slimshady$candidateall_diff <- slimshady$candidate_percent-slimshady$all_percent
slimshady <- slimshady[order(slimshady$candidateall_diff),]


slimshady$odds_ratio <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$candidate_count[x],(slimshady$candidate_count[x]/(slimshady$candidate_percent[x]/100))*(1-slimshady$candidate_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$candidate_count[x],(slimshady$candidate_count[x]/(slimshady$candidate_percent[x]/100))*(1-slimshady$candidate_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$candidate_count[x],(slimshady$candidate_count[x]/(slimshady$candidate_percent[x]/100))*(1-slimshady$candidate_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$candidate_count[x],(slimshady$candidate_count[x]/(slimshady$candidate_percent[x]/100))*(1-slimshady$candidate_percent[x]/100),
                             slimshady$all_count[x],(slimshady$all_count[x]/(slimshady$all_percent[x]/100))*(1-slimshady$all_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
slimshady$logodds <- log(slimshady$odds_ratio)
slimshady <- slimshady[-which(slimshady$go_term=="biological_process"),]
slimshady <- slimshady[-which(slimshady$pval>0.05),]

rm(myCollection,slim,fl,myIds)

#fixing abbreviated GO names
slimshady$go_term<-as.character(slimshady$go_term)
slimshady$go_terms[5]<-"nucleobase-containing compound catabolic process"
slimshady$go_terms[15]<-"generation of precursor metabolites and energy"
slimshady$go_terms[18]<-"cellular nitrogen compound metabolic process"
slimshady$go_terms[19]<-"cellular amino acid metabolic process"
slimshady$go_terms[33]<-"anatomical structure formation involved in morphogenesis"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "sandybrown",mid="white",high = "firebrick4",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("depleted in human lethal genes","intermediate","enriched in human lethal genes"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/go_heatmap.pdf",height=18, width=18, units='cm')


#doing bar graph instead
ggplot(slimshady) +
  geom_bar( aes(x=go_term, y=odds_ratio,fill=odds_ratio ), stat="identity") +
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10")+scale_fill_gradient2(low = "sandybrown", mid = "white",
                                                         high = "firebrick4",trans="log")+  
  theme(axis.text.x = element_text(angle=60, hjust=1))+coord_flip()+
  ylab("Odds Ratio \n (Cell essential OMIM vs Mouse lethal OMIM)")
ggsave("Analysis/Sandra_Figures/Figs/go_bar.pdf",height=18, width=18, units='cm')

