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


slimshady <- slimshady[-which(slimshady$go_term=="biological_process"),]

#only keep 4 most significant terms for dominant, and 4 most significant terms for recessive
slimshady<-slimshady[which((slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio>1),]$pval),3))|
                             (slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio<1),]$pval),3))),]
slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
#fixing abbreviated GO names
slimshady$go_term<-as.character(slimshady$go_term)
slimshady$go_term[which(slimshady$go_term=="anatomical structure formation invo...")]<-"anatomical structure morphogenesis"
slimshady$go_term[which(slimshady$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_BP_poster.pdf",height=18, width=18, units='cm')



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



slimshady <- slimshady[-which(slimshady$go_term=="molecular_function"),]

#only keep 4 most significant terms for dominant, and 4 most significant terms for recessive
slimshady<-slimshady[which((slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio>1),]$pval),3))|
                             (slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio<1),]$pval),3))),]
slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
#fixing abbreviated GO names
slimshady$go_term<-as.character(slimshady$go_term)
slimshady$go_term[which(slimshady$go_term=="DNA binding transcription factor ac...")]<-"DNA binding transcription factor"
slimshady$go_term[which(slimshady$go_term=="transferase activity, transferring ...")]<-" transferase activity"
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)

#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_MF_poster.pdf",height=18, width=18, units='cm')






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


slimshady <- slimshady[-which(slimshady$go_term=="cellular_component"),]

#only keep 4 most significant terms for dominant, and 4 most significant terms for recessive
slimshady<-slimshady[which((slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio>1),]$pval),3))|
                             (slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio<1),]$pval),3))),]
slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))

slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/human_lethal_domrec_CC_poster.pdf",height=18, width=18, units='cm')








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




#only keep 4 most significant terms for dominant, and 4 most significant terms for recessive
slimshady<-slimshady[which((slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio>1),]$pval),3))|
                             (slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio<1),]$pval),3))),]
slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))

#fixing abbreviated GO names
slimshady$go_term<-as.character(slimshady$go_term)
slimshady$go_term[which(slimshady$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshady$go_term[which(slimshady$go_term=="cellular protein modification proce...")]<-"cellular protein modification process"

slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("Enriched in PUTATIVE RECESSIVE human lethal genes","","enriched in PUTATIVE DOMINANT human lethal genes"))+  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/PUTATIVE_domrec_BP_poster.pdf",height=18, width=18, units='cm')








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



slimshady <- slimshady[-which(slimshady$go_term=="molecular_function"),]

#only keep 4 most significant terms for dominant, and 4 most significant terms for recessive
slimshady<-slimshady[which((slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio>1),]$pval),3))|
                             (slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio<1),]$pval),3))),]
slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))
slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)

#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",
                       breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/PUTATIVE_domrec_MF_poster_legend.pdf",height=18, width=18, units='cm')







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


slimshady <- slimshady[-which(slimshady$go_term=="cellular_component"),]

#only keep 4 most significant terms for dominant, and 4 most significant terms for recessive
slimshady<-slimshady[which((slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio>1),]$pval),3))|
                             (slimshady$pval%in%head(sort(slimshady[which(slimshady$odds_ratio<1),]$pval),3))),]
slimshady <- slimshady[order(slimshady$odds_ratio),]

slimshady$gene <- rep("gene",length(slimshady$go_term))

slimshady$go_term <- factor(slimshady$go_term, levels = slimshady$go_term)
#plotting heat map
ggplot(data = slimshady, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +scale_fill_gradient2(low = "#2c5daa",mid="white",high = "#e32026",trans="log",breaks=c(min(slimshady$odds_ratio),1,max(slimshady$odds_ratio)),labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=20))+scale_y_discrete(expand = c(0, 0))+scale_x_discrete(expand = c(0, 0))
ggsave("Analysis/Sandra_Figures/Figs/PUTATIVE_domrec_CC_poster.pdf",height=18, width=18, units='cm')







