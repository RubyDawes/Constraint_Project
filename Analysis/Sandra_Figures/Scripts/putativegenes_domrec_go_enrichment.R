library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############PUTATIVE lethal genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshady$rec_count <- a$Count
slimshady$rec_percent <- a$Percent
rm(a)

slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$dom_percent[which(slimshady$dom_count==0)]<-slimshady$dom_percent[which(slimshady$dom_count==1)][1]
slimshady$rec_percent[which(slimshady$rec_percent==0)]<-slimshady$rec_percent[which(slimshady$rec_count==1)][1]
slimshady$dom_count[which(slimshady$dom_count==0)]<-1
slimshady$rec_count[which(slimshady$rec_count==0)]<-1

slimshady$domrec_diff <- slimshady$dom_percent-slimshady$rec_percent


slimshady$odds_ratio <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
slimshady$logodds <- log(slimshady$odds_ratio)
slimshady <- slimshady[-which(slimshady$go_term=="biological_process"),]
slimshady <- slimshady[-which(slimshady$pval>0.005),]

#fixing abbreviated GO names
slimshady$go_term<-as.character(slimshady$go_term)
slimshady$go_term[which(slimshady$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshady$go_term[which(slimshady$go_term=="cellular protein modification proce...")]<-"cellular protein modification process"
slimshady$go_term[which(slimshady$go_term=="cellular nitrogen compound metaboli...")]<-"cellular nitrogen compound metabolic process"
slimshady$go_term[which(slimshady$go_term=="generation of precursor metabolites...")]<-"generation of precursor metabolites and energy"

slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "#beaed4",mid="white",high = "#7fc97f",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("Enriched in PUTATIVE RECESSIVE human lethal genes","","enriched in PUTATIVE DOMINANT human lethal genes"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/PUTATIVE_domrec_BP_go_heatmap_0.005.pdf",height=18, width=35, units='cm')






############PUTATIVE lethal genes dominant vs recessive MF################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "MF")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "MF")

slimshady$rec_count <- a$Count
slimshady$rec_percent <- a$Percent
rm(a)

slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$dom_percent[which(slimshady$dom_count==0)]<-slimshady$dom_percent[which(slimshady$dom_count==1)][1]
slimshady$rec_percent[which(slimshady$rec_percent==0)]<-slimshady$rec_percent[which(slimshady$rec_count==2)][1]/2
slimshady$dom_count[which(slimshady$dom_count==0)]<-1
slimshady$rec_count[which(slimshady$rec_count==0)]<-1

slimshady$domrec_diff <- slimshady$dom_percent-slimshady$rec_percent


slimshady$odds_ratio <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
slimshady$logodds <- log(slimshady$odds_ratio)
slimshady <- slimshady[-which(slimshady$go_term=="molecular_function"),]
slimshady <- slimshady[-which(slimshady$pval>0.05),]
#fixing abbreviated GO names
slimshady$go_term<-as.character(slimshady$go_term)

slimshady$go_term[which(slimshady$go_term=="transcription factor activity, prot...")]<-"transcription factor activity, protein binding"
slimshady$go_term[which(slimshady$go_term=="hydrolase activity, acting on carbo...")]<-"hydrolase activity"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)

#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "#beaed4",mid="white",high = "#7fc97f",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("Enriched in PUTATIVE RECESSIVE human lethal genes","","enriched in PUTATIVE DOMINANT human lethal genes"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/PUTATIVE_domrec_MF_go_heatmap.pdf",height=18, width=30, units='cm')







############PUTATIVE lethal genes dominant vs recessive CC################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "CC")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "CC")

slimshady$rec_count <- a$Count
slimshady$rec_percent <- a$Percent
rm(a)

slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$dom_percent[which(slimshady$dom_count==0)]<-slimshady$dom_percent[which(slimshady$dom_count==2)][1]/2
slimshady$rec_percent[which(slimshady$rec_percent==0)]<-slimshady$rec_percent[which(slimshady$rec_count==2)][1]/2
slimshady$dom_count[which(slimshady$dom_count==0)]<-1
slimshady$rec_count[which(slimshady$rec_count==0)]<-1

slimshady$domrec_diff <- slimshady$dom_percent-slimshady$rec_percent


slimshady$odds_ratio <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
slimshady$logodds <- log(slimshady$odds_ratio)
slimshady <- slimshady[-which(slimshady$go_term=="cellular_component"),]
slimshady <- slimshady[-which(slimshady$pval>0.05),]
#fixing abbreviated GO names
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "#beaed4",mid="white",high = "#7fc97f",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("Enriched in PUTATIVE RECESSIVE human lethal genes","","enriched in PUTATIVE DOMINANT human lethal genes"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/PUTATIVE_domrec_CC_go_heatmap.pdf",height=18, width=30, units='cm')




