library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############lethal genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshady$rec_count <- a$Count
slimshady$rec_percent <- a$Percent
rm(a)

#all genes go slim annotations- human lethal INBETWEENERS
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AR,AD"|universe_df$Inheritance_pattern=="MT,AR,AD"|universe_df$Inheritance_pattern=="MT,XLd,AR"|universe_df$Inheritance_pattern=="XLd,XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshady$inb_count <- a$Count
slimshady$inb_percent <- a$Percent
rm(a)

slimshady<-slimshady[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshady$dom_percent[which(slimshady$dom_count==0)]<-slimshady$dom_percent[which(slimshady$dom_count==1)][1]
slimshady$rec_percent[which(slimshady$rec_percent==0)]<-slimshady$rec_percent[which(slimshady$rec_count==2)][1]/2
slimshady$inb_percent[which(slimshady$inb_percent==0)]<-slimshady$rec_percent[which(slimshady$inb_count==1)][1]
slimshady$dom_count[which(slimshady$dom_count==0)]<-1
slimshady$rec_count[which(slimshady$rec_count==0)]<-1
slimshady$inb_count[which(slimshady$inb_count==0)]<-1

slimshady$inbrec_diff <- slimshady$rec_percent-slimshady$inb_percent
slimshady$inbdom_diff <- slimshady$dom_percent-slimshady$inb_percent


slimshady$odds_ratio_inbrec <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower_inbrec <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher_inbrec <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval_inbrec <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100),
                             slimshady$rec_count[x],(slimshady$rec_count[x]/(slimshady$rec_percent[x]/100))*(1-slimshady$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

slimshady$odds_ratio_inbdom <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshady$confint_lower_inbdom <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshady$confint_higher_inbdom <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshady$pval_inbdom <- unlist(lapply(1:length(slimshady$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshady$dom_count[x],(slimshady$dom_count[x]/(slimshady$dom_percent[x]/100))*(1-slimshady$dom_percent[x]/100),
                             slimshady$inb_count[x],(slimshady$inb_count[x]/(slimshady$inb_percent[x]/100))*(1-slimshady$inb_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))

#INBETWEENERS VS DOMINANT
slimshady <- slimshady[order(slimshady$odds_ratio_inbdom),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
slimshady <- slimshady[-which(slimshady$go_term=="biological_process"),]

#fixing abbreviated GO names
#slimshady$go_term<-as.character(slimshady$go_term)

#slimshady$go_term[which(slimshady$go_term=="nucleobase-containing compound cata...")]<-"nucleobase-containing compound catabolic process"
#slimshady$go_term[which(slimshady$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
inbdom<- slimshady[-which(slimshady$pval_inbdom>0.05),]
ggplot(data = inbdom, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_inbdom,width=1)) +
  scale_fill_gradient2(low = "pink",mid="white",high = "green",trans="log",
                       breaks=c(min(slimshady$odds_ratio_inbdom),1,max(slimshady$odds_ratio_inbdom)),
                       labels=c("Enriched in INBETWEENER human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_inbdom_BP_go_heatmap.pdf",height=18, width=35, units='cm')

#INBETWEENERS VS RECESSIVE
slimshady <- slimshady[order(slimshady$odds_ratio_inbrec),]
#fixing abbreviated GO names
#slimshady$go_term<-as.character(slimshady$go_term)

#slimshady$go_term[which(slimshady$go_term=="nucleobase-containing compound cata...")]<-"nucleobase-containing compound catabolic process"
#slimshady$go_term[which(slimshady$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
inbrec<- slimshady[-which(slimshady$pval_inbrec>0.05),]
ggplot(data = inbrec, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_inbrec,width=1)) +
  scale_fill_gradient2(low = "pink",mid="white",high = "green",trans="log",
                       breaks=c(min(slimshady$odds_ratio_inbrec),1,max(slimshady$odds_ratio_inbrec)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in INBETWEENER human lethal genes"))+
  bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_inbdom_BP_go_heatmap.pdf",height=18, width=35, units='cm')






############lethal genes dominant vs recessive MF################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "MF")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "MF")

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
slimshady <- slimshady[-which(slimshady$go_term=="molecular_function"),]
slimshady <- slimshady[-which(slimshady$pval>0.05),]
#fixing abbreviated GO names
slimshady$go_term<-as.character(slimshady$go_term)

slimshady$go_term[which(slimshady$go_term=="DNA binding transcription factor ac...")]<-"DNA-binding transcription factor activity"
slimshady$go_term[which(slimshady$go_term=="transferase activity, transferring ...")]<-"transferase"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)

#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "pink",mid="white",high = "green",trans="log",
                       breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_MF_go_heatmap.pdf",height=18, width=30, units='cm')






############lethal genes dominant vs recessive CC################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "CC")

slimshady <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "CC")

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
slimshady <- slimshady[-which(slimshady$go_term=="cellular_component"),]
slimshady <- slimshady[-which(slimshady$pval>0.05),]
#fixing abbreviated GO names
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "pink",mid="white",high = "green",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+bar_theme()
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_CC_go_heatmap.pdf",height=18, width=30, units='cm')




