library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############mouse lethal vs mouse non-lethal################
#all genes go slim annotations- mouse lethal
all <- unlist(universe_df$go_terms[which(universe_df$lethal_mouse=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,lethal_count = aa$Count,lethal_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- mouse non-lethal
myIds <- unlist(universe_df$go_terms[which(universe_df$lethal_mouse=="N"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshady$nonlethal_count <- a$Count
slimshady$nonlethal_percent <- a$Percent
rm(a)


slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$nonlethal_percent[which(slimshady$nonlethal_percent==0)]<-slimshady$nonlethal_percent[which(slimshady$nonlethal_count==6)][1]/6
slimshady$nonlethal_count[which(slimshady$nonlethal_count==0)]<-1

slimshady$lethalall_diff <- slimshady$lethal_percent-slimshady$nonlethal_percent
slimshady <- slimshady[order(slimshady$lethalall_diff),]


slimshady$odds_ratio <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$nonlethal_count[x],(slimshady$nonlethal_count[x]/(slimshady$nonlethal_percent[x]/100))*(1-slimshady$nonlethal_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$nonlethal_count[x],(slimshady$nonlethal_count[x]/(slimshady$nonlethal_percent[x]/100))*(1-slimshady$nonlethal_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$nonlethal_count[x],(slimshady$nonlethal_count[x]/(slimshady$nonlethal_percent[x]/100))*(1-slimshady$nonlethal_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$lethal_count[x],(slimshady$lethal_count[x]/(slimshady$lethal_percent[x]/100))*(1-slimshady$lethal_percent[x]/100),
                             slimshady$nonlethal_count[x],(slimshady$nonlethal_count[x]/(slimshady$nonlethal_percent[x]/100))*(1-slimshady$nonlethal_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
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
slimshady$go_term[which(slimshady$go_term=="nucleobase-containing compound cata...")]<-"nucleobase-containing compound catabolic process"
slimshady$go_term[which(slimshady$go_term=="cellular nitrogen compound metaboli...")]<-"cellular nitrogen compound metabolic process"
slimshady$go_term[which(slimshady$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshady$go_term[which(slimshady$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshady$go_term[which(slimshady$go_term=="cytoskeleton-dependent intracellula...")]<-"cytoskeleton-dependent intracellular transport"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "steelblue3",mid="white",high = "black",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("enriched in mouse non-lethal genes","","enriched in mouse lethal genes"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/mouselethal_go_heatmap.pdf",height=18, width=22, units='cm')

############pie of mouse lethal genes##############
mouse_lethal_pie <- data.frame(group=c("Lethal     ","Non-Lethal     "),
                               value=c(length(which(universe_df$lethal_mouse=="Y")),
                                       length(which(universe_df$lethal_mouse=="N"))))
mouse_lethal_pie$group <- factor(mouse_lethal_pie$group, levels = mouse_lethal_pie$group)
c <- ggplot(mouse_lethal_pie, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("black","steelblue3"))
c<- c+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
ggsave("Analysis/Sandra_Figures/Figs/mouselethal_pie.pdf",height=18, width=18, units='cm')

