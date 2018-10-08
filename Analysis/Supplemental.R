source("plot_functions.R")

###Figure 5: Ontology analysis of known prenatal/infantile lethal genes and putative prenatal/ infantile lethal genes#####
library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############lethal genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshadyBP <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshadyBP$rec_count <- a$Count
slimshadyBP$rec_percent <- a$Percent
rm(a)

slimshadyBP<-slimshadyBP[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyBP$dom_percent[which(slimshadyBP$dom_count==0)]<-slimshadyBP$dom_percent[which(slimshadyBP$dom_count==1)][1]
slimshadyBP$dom_count[which(slimshadyBP$dom_count==0)]<-1

slimshadyBP$domrec_diff <- slimshadyBP$dom_percent-slimshadyBP$rec_percent


slimshadyBP$odds_ratio <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyBP$confint_lower <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyBP$confint_higher <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyBP$pval <- unlist(lapply(1:length(slimshadyBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBP$dom_count[x],(slimshadyBP$dom_count[x]/(slimshadyBP$dom_percent[x]/100))*(1-slimshadyBP$dom_percent[x]/100),
                             slimshadyBP$rec_count[x],(slimshadyBP$rec_count[x]/(slimshadyBP$rec_percent[x]/100))*(1-slimshadyBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyBP <- slimshadyBP[-which(slimshadyBP$go_term=="biological_process"),]

slimshadyBP <- slimshadyBP[order(slimshadyBP$odds_ratio),]

slimshadyBP$gene <- rep("gene",length(slimshadyBP$go_term))
#fixing abbreviated GO names
slimshadyBP$go_term<-as.character(slimshadyBP$go_term)
slimshadyBP$go_term[which(slimshadyBP$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"


slimshadyBP$go_term <- factor(slimshadyBP$go_term, levels = slimshadyBP$go_term)

# taking 10 most significantly different terms
slimshadyBP <- slimshadyBP[order(slimshadyBP$pval),]
slimshadyBP <- slimshadyBP[1:8,]
slimshadyBP <- slimshadyBP[order(slimshadyBP$odds_ratio),]
#plotting heat map
ap<-ggplot(data = slimshadyBP, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45" ,trans="log",
                       breaks=c(min(slimshadyBP$odds_ratio),1,max(slimshadyBP$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 20, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

############lethal genes dominant vs recessive CC################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "CC")

slimshadyCC <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "CC")

slimshadyCC$rec_count <- a$Count
slimshadyCC$rec_percent <- a$Percent
rm(a)

slimshadyCC<-slimshadyCC[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyCC$dom_percent[which(slimshadyCC$dom_count==0)]<-slimshadyCC$dom_percent[which(slimshadyCC$dom_count==1)][1]
slimshadyCC$rec_percent[which(slimshadyCC$rec_percent==0)]<-slimshadyCC$rec_percent[which(slimshadyCC$rec_count==6)][1]/6
slimshadyCC$dom_count[which(slimshadyCC$dom_count==0)]<-1
slimshadyCC$rec_count[which(slimshadyCC$rec_count==0)]<-1

slimshadyCC$domrec_diff <- slimshadyCC$dom_percent-slimshadyCC$rec_percent


slimshadyCC$odds_ratio <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyCC$confint_lower <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyCC$confint_higher <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyCC$pval <- unlist(lapply(1:length(slimshadyCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCC$dom_count[x],(slimshadyCC$dom_count[x]/(slimshadyCC$dom_percent[x]/100))*(1-slimshadyCC$dom_percent[x]/100),
                             slimshadyCC$rec_count[x],(slimshadyCC$rec_count[x]/(slimshadyCC$rec_percent[x]/100))*(1-slimshadyCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyCC <- slimshadyCC[-which(slimshadyCC$go_term=="cellular_component"),]

slimshadyCC <- slimshadyCC[order(slimshadyCC$odds_ratio),]

slimshadyCC$gene <- rep("gene",length(slimshadyCC$go_term))

slimshadyCC$go_term <- factor(slimshadyCC$go_term, levels = slimshadyCC$go_term)

# taking 10 most significantly different terms
slimshadyCC <- slimshadyCC[order(slimshadyCC$pval),]
slimshadyCC <- slimshadyCC[1:6,]
slimshadyCC <- slimshadyCC[order(slimshadyCC$odds_ratio),]

#plotting heat map
b<-ggplot(data = slimshadyCC, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45" ,trans="log",
                       breaks=c(min(slimshadyCC$odds_ratio),1,max(slimshadyCC$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))








###########lethal genes dominant vs recessive MF################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "MF")

slimshadyMF <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_B=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "MF")

slimshadyMF$rec_count <- a$Count
slimshadyMF$rec_percent <- a$Percent
rm(a)

slimshadyMF<-slimshadyMF[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyMF$dom_percent[which(slimshadyMF$dom_count==0)]<-slimshadyMF$dom_percent[which(slimshadyMF$dom_count==1)][1]
slimshadyMF$rec_percent[which(slimshadyMF$rec_percent==0)]<-slimshadyMF$rec_percent[which(slimshadyMF$rec_count==1)][1]
slimshadyMF$dom_count[which(slimshadyMF$dom_count==0)]<-1
slimshadyMF$rec_count[which(slimshadyMF$rec_count==0)]<-1

slimshadyMF$domrec_diff <- slimshadyMF$dom_percent-slimshadyMF$rec_percent


slimshadyMF$odds_ratio <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyMF$confint_lower <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyMF$confint_higher <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyMF$pval <- unlist(lapply(1:length(slimshadyMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMF$dom_count[x],(slimshadyMF$dom_count[x]/(slimshadyMF$dom_percent[x]/100))*(1-slimshadyMF$dom_percent[x]/100),
                             slimshadyMF$rec_count[x],(slimshadyMF$rec_count[x]/(slimshadyMF$rec_percent[x]/100))*(1-slimshadyMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyMF <- slimshadyMF[-which(slimshadyMF$go_term=="molecular_function"),]

slimshadyMF <- slimshadyMF[order(slimshadyMF$odds_ratio),]

slimshadyMF$gene <- rep("gene",length(slimshadyMF$go_term))


# taking 10 most significantly different terms
slimshadyMF <- slimshadyMF[order(slimshadyMF$pval),]
slimshadyMF <- slimshadyMF[1:8,]
slimshadyMF <- slimshadyMF[order(slimshadyMF$odds_ratio),]

#fixing abbreviated GO names
slimshadyMF$go_term<-as.character(slimshadyMF$go_term)
slimshadyMF$go_term[which(slimshadyMF$go_term=="DNA binding transcription factor ac...")]<-"DNA binding transcription factor"
slimshadyMF$go_term[which(slimshadyMF$go_id=="GO:0016757")]<-"glycosyltransferase activity"



slimshadyMF$go_term <- factor(slimshadyMF$go_term, levels = slimshadyMF$go_term)

#plotting heat map
c<-ggplot(data = slimshadyMF, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45" ,trans="log",
                       breaks=c(min(slimshadyMF$odds_ratio),1,max(slimshadyMF$odds_ratio)-0.2),
                       labels=c("RECESSIVE           ","","DOMINANT        "))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 3, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))




#saving plot####
ggarrange(ap,b,c,widths=c(1,1.1,1.65),ncol=3)
ggsave("output/Figures/Supp.pdf",height=15, width=26, units='cm')


rm(ap,b,c,d,e,f,slim,myCollection,slimshadyBP,slimshadyBPput,slimshadyCC,slimshadyCCput,slimshadyMF,slimshadyMFput,fl,myIds)