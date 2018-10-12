source("plot_functions.R")

###Figure 5: Ontology analysis of known prenatal/infantile lethal genes and putative prenatal/ infantile lethal genes#####
library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############omim genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshadyOMIMOMIMBP <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshadyOMIMOMIMBP$rec_count <- a$Count
slimshadyOMIMOMIMBP$rec_percent <- a$Percent
rm(a)

slimshadyOMIMOMIMBP<-slimshadyOMIMOMIMBP[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyOMIMOMIMBP$dom_percent[which(slimshadyOMIMOMIMBP$dom_count==0)]<-slimshadyOMIMOMIMBP$dom_percent[which(slimshadyOMIMOMIMBP$dom_count==1)][1]
slimshadyOMIMOMIMBP$dom_count[which(slimshadyOMIMOMIMBP$dom_count==0)]<-1

slimshadyOMIMOMIMBP$domrec_diff <- slimshadyOMIMOMIMBP$dom_percent-slimshadyOMIMOMIMBP$rec_percent


slimshadyOMIMOMIMBP$odds_ratio <- unlist(lapply(1:length(slimshadyOMIMOMIMBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMOMIMBP$dom_count[x],(slimshadyOMIMOMIMBP$dom_count[x]/(slimshadyOMIMBP$dom_percent[x]/100))*(1-slimshadyOMIMBP$dom_percent[x]/100),
                             slimshadyOMIMBP$rec_count[x],(slimshadyOMIMBP$rec_count[x]/(slimshadyOMIMBP$rec_percent[x]/100))*(1-slimshadyOMIMBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyOMIMBP$confint_lower <- unlist(lapply(1:length(slimshadyOMIMBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMBP$dom_count[x],(slimshadyOMIMBP$dom_count[x]/(slimshadyOMIMBP$dom_percent[x]/100))*(1-slimshadyOMIMBP$dom_percent[x]/100),
                             slimshadyOMIMBP$rec_count[x],(slimshadyOMIMBP$rec_count[x]/(slimshadyOMIMBP$rec_percent[x]/100))*(1-slimshadyOMIMBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyOMIMBP$confint_higher <- unlist(lapply(1:length(slimshadyOMIMBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMBP$dom_count[x],(slimshadyOMIMBP$dom_count[x]/(slimshadyOMIMBP$dom_percent[x]/100))*(1-slimshadyOMIMBP$dom_percent[x]/100),
                             slimshadyOMIMBP$rec_count[x],(slimshadyOMIMBP$rec_count[x]/(slimshadyOMIMBP$rec_percent[x]/100))*(1-slimshadyOMIMBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyOMIMBP$pval <- unlist(lapply(1:length(slimshadyOMIMBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMBP$dom_count[x],(slimshadyOMIMBP$dom_count[x]/(slimshadyOMIMBP$dom_percent[x]/100))*(1-slimshadyOMIMBP$dom_percent[x]/100),
                             slimshadyOMIMBP$rec_count[x],(slimshadyOMIMBP$rec_count[x]/(slimshadyOMIMBP$rec_percent[x]/100))*(1-slimshadyOMIMBP$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyOMIMBP <- slimshadyOMIMBP[-which(slimshadyOMIMBP$go_term=="biological_process"),]

slimshadyOMIMBP <- slimshadyOMIMBP[order(slimshadyOMIMBP$odds_ratio),]

slimshadyOMIMBP$gene <- rep("gene",length(slimshadyOMIMBP$go_term))
#fixing abbreviated GO names
slimshadyOMIMBP$go_term<-as.character(slimshadyOMIMBP$go_term)
slimshadyOMIMBP$go_term[which(slimshadyOMIMBP$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"


slimshadyOMIMBP$go_term <- factor(slimshadyOMIMBP$go_term, levels = slimshadyOMIMBP$go_term)

# taking 10 most significantly different terms
slimshadyOMIMBP <- slimshadyOMIMBP[order(slimshadyOMIMBP$pval),]
slimshadyOMIMBP <- slimshadyOMIMBP[1:8,]
slimshadyOMIMBP <- slimshadyOMIMBP[order(slimshadyOMIMBP$odds_ratio),]
#plotting heat map
api<-ggplot(data = slimshadyOMIMBP, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45" ,trans="log",
                       breaks=c(min(slimshadyOMIMBP$odds_ratio),1,max(slimshadyOMIMBP$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 20, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

############omim genes dominant vs recessive CC################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "CC")

slimshadyOMIMCC <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "CC")

slimshadyOMIMCC$rec_count <- a$Count
slimshadyOMIMCC$rec_percent <- a$Percent
rm(a)

slimshadyOMIMCC<-slimshadyOMIMCC[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyOMIMCC$dom_percent[which(slimshadyOMIMCC$dom_count==0)]<-slimshadyOMIMCC$dom_percent[which(slimshadyOMIMCC$dom_count==1)][1]
slimshadyOMIMCC$rec_percent[which(slimshadyOMIMCC$rec_percent==0)]<-slimshadyOMIMCC$rec_percent[which(slimshadyOMIMCC$rec_count==6)][1]/6
slimshadyOMIMCC$dom_count[which(slimshadyOMIMCC$dom_count==0)]<-1
slimshadyOMIMCC$rec_count[which(slimshadyOMIMCC$rec_count==0)]<-1

slimshadyOMIMCC$domrec_diff <- slimshadyOMIMCC$dom_percent-slimshadyOMIMCC$rec_percent


slimshadyOMIMCC$odds_ratio <- unlist(lapply(1:length(slimshadyOMIMCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMCC$dom_count[x],(slimshadyOMIMCC$dom_count[x]/(slimshadyOMIMCC$dom_percent[x]/100))*(1-slimshadyOMIMCC$dom_percent[x]/100),
                             slimshadyOMIMCC$rec_count[x],(slimshadyOMIMCC$rec_count[x]/(slimshadyOMIMCC$rec_percent[x]/100))*(1-slimshadyOMIMCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyOMIMCC$confint_lower <- unlist(lapply(1:length(slimshadyOMIMCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMCC$dom_count[x],(slimshadyOMIMCC$dom_count[x]/(slimshadyOMIMCC$dom_percent[x]/100))*(1-slimshadyOMIMCC$dom_percent[x]/100),
                             slimshadyOMIMCC$rec_count[x],(slimshadyOMIMCC$rec_count[x]/(slimshadyOMIMCC$rec_percent[x]/100))*(1-slimshadyOMIMCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyOMIMCC$confint_higher <- unlist(lapply(1:length(slimshadyOMIMCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMCC$dom_count[x],(slimshadyOMIMCC$dom_count[x]/(slimshadyOMIMCC$dom_percent[x]/100))*(1-slimshadyOMIMCC$dom_percent[x]/100),
                             slimshadyOMIMCC$rec_count[x],(slimshadyOMIMCC$rec_count[x]/(slimshadyOMIMCC$rec_percent[x]/100))*(1-slimshadyOMIMCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyOMIMCC$pval <- unlist(lapply(1:length(slimshadyOMIMCC$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMCC$dom_count[x],(slimshadyOMIMCC$dom_count[x]/(slimshadyOMIMCC$dom_percent[x]/100))*(1-slimshadyOMIMCC$dom_percent[x]/100),
                             slimshadyOMIMCC$rec_count[x],(slimshadyOMIMCC$rec_count[x]/(slimshadyOMIMCC$rec_percent[x]/100))*(1-slimshadyOMIMCC$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyOMIMCC <- slimshadyOMIMCC[-which(slimshadyOMIMCC$go_term=="cellular_component"),]

slimshadyOMIMCC <- slimshadyOMIMCC[order(slimshadyOMIMCC$odds_ratio),]

slimshadyOMIMCC$gene <- rep("gene",length(slimshadyOMIMCC$go_term))

slimshadyOMIMCC$go_term <- factor(slimshadyOMIMCC$go_term, levels = slimshadyOMIMCC$go_term)

# taking 10 most significantly different terms
slimshadyOMIMCC <- slimshadyOMIMCC[order(slimshadyOMIMCC$pval),]
slimshadyOMIMCC <- slimshadyOMIMCC[1:6,]
slimshadyOMIMCC <- slimshadyOMIMCC[order(slimshadyOMIMCC$odds_ratio),]

#plotting heat map
bi<-ggplot(data = slimshadyOMIMCC, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45" ,trans="log",
                       breaks=c(min(slimshadyOMIMCC$odds_ratio),1,max(slimshadyOMIMCC$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))








###########omim genes dominant vs recessive MF################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "MF")

slimshadyOMIMMF <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_B=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "MF")

slimshadyOMIMMF$rec_count <- a$Count
slimshadyOMIMMF$rec_percent <- a$Percent
rm(a)

slimshadyOMIMMF<-slimshadyOMIMMF[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyOMIMMF$dom_percent[which(slimshadyOMIMMF$dom_count==0)]<-slimshadyOMIMMF$dom_percent[which(slimshadyOMIMMF$dom_count==1)][1]
slimshadyOMIMMF$rec_percent[which(slimshadyOMIMMF$rec_percent==0)]<-slimshadyOMIMMF$rec_percent[which(slimshadyOMIMMF$rec_count==1)][1]
slimshadyOMIMMF$dom_count[which(slimshadyOMIMMF$dom_count==0)]<-1
slimshadyOMIMMF$rec_count[which(slimshadyOMIMMF$rec_count==0)]<-1

slimshadyOMIMMF$domrec_diff <- slimshadyOMIMMF$dom_percent-slimshadyOMIMMF$rec_percent


slimshadyOMIMMF$odds_ratio <- unlist(lapply(1:length(slimshadyOMIMMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMMF$dom_count[x],(slimshadyOMIMMF$dom_count[x]/(slimshadyOMIMMF$dom_percent[x]/100))*(1-slimshadyOMIMMF$dom_percent[x]/100),
                             slimshadyOMIMMF$rec_count[x],(slimshadyOMIMMF$rec_count[x]/(slimshadyOMIMMF$rec_percent[x]/100))*(1-slimshadyOMIMMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyOMIMMF$confint_lower <- unlist(lapply(1:length(slimshadyOMIMMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMMF$dom_count[x],(slimshadyOMIMMF$dom_count[x]/(slimshadyOMIMMF$dom_percent[x]/100))*(1-slimshadyOMIMMF$dom_percent[x]/100),
                             slimshadyOMIMMF$rec_count[x],(slimshadyOMIMMF$rec_count[x]/(slimshadyOMIMMF$rec_percent[x]/100))*(1-slimshadyOMIMMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyOMIMMF$confint_higher <- unlist(lapply(1:length(slimshadyOMIMMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMMF$dom_count[x],(slimshadyOMIMMF$dom_count[x]/(slimshadyOMIMMF$dom_percent[x]/100))*(1-slimshadyOMIMMF$dom_percent[x]/100),
                             slimshadyOMIMMF$rec_count[x],(slimshadyOMIMMF$rec_count[x]/(slimshadyOMIMMF$rec_percent[x]/100))*(1-slimshadyOMIMMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyOMIMMF$pval <- unlist(lapply(1:length(slimshadyOMIMMF$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMMF$dom_count[x],(slimshadyOMIMMF$dom_count[x]/(slimshadyOMIMMF$dom_percent[x]/100))*(1-slimshadyOMIMMF$dom_percent[x]/100),
                             slimshadyOMIMMF$rec_count[x],(slimshadyOMIMMF$rec_count[x]/(slimshadyOMIMMF$rec_percent[x]/100))*(1-slimshadyOMIMMF$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyOMIMMF <- slimshadyOMIMMF[-which(slimshadyOMIMMF$go_term=="molecular_function"),]

slimshadyOMIMMF <- slimshadyOMIMMF[order(slimshadyOMIMMF$odds_ratio),]

slimshadyOMIMMF$gene <- rep("gene",length(slimshadyOMIMMF$go_term))


# taking 10 most significantly different terms
slimshadyOMIMMF <- slimshadyOMIMMF[order(slimshadyOMIMMF$pval),]
slimshadyOMIMMF <- slimshadyOMIMMF[1:8,]
slimshadyOMIMMF <- slimshadyOMIMMF[order(slimshadyOMIMMF$odds_ratio),]

#fixing abbreviated GO names
slimshadyOMIMMF$go_term<-as.character(slimshadyOMIMMF$go_term)
slimshadyOMIMMF$go_term[which(slimshadyOMIMMF$go_term=="DNA binding transcription factor ac...")]<-"DNA binding transcription factor"
slimshadyOMIMMF$go_term[which(slimshadyOMIMMF$go_id=="GO:0016757")]<-"glycosyltransferase activity"



slimshadyOMIMMF$go_term <- factor(slimshadyOMIMMF$go_term, levels = slimshadyOMIMMF$go_term)

#plotting heat map
ci<-ggplot(data = slimshadyOMIMMF, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45" ,trans="log",
                       breaks=c(min(slimshadyOMIMMF$odds_ratio),1,max(slimshadyOMIMMF$odds_ratio)-0.2),
                       labels=c("RECESSIVE           ","","DOMINANT        "))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 3, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))




#saving plot####
ggarrange(ap,b,c,widths=c(1,1.1,1.65),ncol=3)
ggsave("output/Figures/Supp.pdf",height=15, width=26, units='cm')


rm(ap,b,c,d,e,f,slim,myCollection,slimshadyOMIMBP,slimshadyOMIMBPput,slimshadyOMIMCC,slimshadyOMIMCCput,slimshadyOMIMMF,slimshadyOMIMMFput,fl,myIds)