source("plot_functions.R")

###Figure 5: Ontology analysis of known prenatal/infantile lethal genes and putative prenatal/ infantile lethal genes#####
library(GSEABase)

fl <- "http://www.geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
############lethal genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshadyBP <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshadyBP$rec_count <- a$Count
slimshadyBP$rec_percent <- a$Percent
rm(a)

slimshadyBP<-slimshadyBP[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyBP$dom_percent[which(slimshadyBP$dom_count==0)]<-slimshadyBP$dom_percent[which(slimshadyBP$dom_count==1)][1]
slimshadyBP$rec_percent[which(slimshadyBP$rec_percent==0)]<-slimshadyBP$rec_percent[which(slimshadyBP$rec_count==2)][1]/2
slimshadyBP$dom_count[which(slimshadyBP$dom_count==0)]<-1
slimshadyBP$rec_count[which(slimshadyBP$rec_count==0)]<-1

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
slimshadyBP$go_term[which(slimshadyBP$go_id=="GO:0034655")]<-"nucleobase-containing compound catabolic process"
slimshadyBP$go_term[which(slimshadyBP$go_id=="GO:0030705")]<-"cytoskeleton-dependent intracellular transport"
slimshadyBP$go_term[which(slimshadyBP$go_id=="GO:0006091")]<-"generation of precursor metabolites and energy"


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
  which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "CC")

slimshadyCC <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "CC")

slimshadyCC$rec_count <- a$Count
slimshadyCC$rec_percent <- a$Percent
rm(a)

slimshadyCC<-slimshadyCC[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyCC$dom_percent[which(slimshadyCC$dom_count==0)]<-slimshadyCC$dom_percent[which(slimshadyCC$dom_count==1)][1]
slimshadyCC$rec_percent[which(slimshadyCC$rec_percent==0)]<-slimshadyCC$rec_percent[which(slimshadyCC$rec_count==1)][1]
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
  which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&lengths(universe_df$go_terms)>0)])
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
#fixing abbreviated GO names
slimshadyMF$go_term<-as.character(slimshadyMF$go_term)
slimshadyMF$go_term[which(slimshadyMF$go_term=="DNA binding transcription factor ac...")]<-"DNA binding transcription factor"
slimshadyMF$go_term[which(slimshadyMF$go_id=="GO:0016757")]<-"glycosyltransferase activity"
slimshadyMF$go_term[which(slimshadyMF$go_id=="GO:0016765")]<-"transferase activity, transferring alkyl or aryl (other than methyl) groups"



slimshadyMF$go_term <- factor(slimshadyMF$go_term, levels = slimshadyMF$go_term)


# taking 10 most significantly different terms
slimshadyMF <- slimshadyMF[order(slimshadyMF$pval),]
slimshadyMF <- slimshadyMF[1:8,]
slimshadyMF <- slimshadyMF[order(slimshadyMF$odds_ratio),]
#plotting heat map
c<-ggplot(data = slimshadyMF, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45" ,trans="log",
                       breaks=c(min(slimshadyMF$odds_ratio),1,max(slimshadyMF$odds_ratio)-0.2),
                       labels=c("RECESSIVE           ","","DOMINANT        "))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 3, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))



############PUTATIVE lethal genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y")&universe_df$constrained=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshadyBPput <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y")&universe_df$constrained=="N"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshadyBPput$rec_count <- a$Count
slimshadyBPput$rec_percent <- a$Percent
rm(a)

slimshadyBPput<-slimshadyBPput[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyBPput$dom_percent[which(slimshadyBPput$dom_count==0)]<-slimshadyBPput$dom_percent[which(slimshadyBPput$dom_count==1)][1]
slimshadyBPput$rec_percent[which(slimshadyBPput$rec_percent==0)]<-slimshadyBPput$rec_percent[which(slimshadyBPput$rec_count==1)][1]
slimshadyBPput$dom_count[which(slimshadyBPput$dom_count==0)]<-1
slimshadyBPput$rec_count[which(slimshadyBPput$rec_count==0)]<-1

slimshadyBPput$domrec_diff <- slimshadyBPput$dom_percent-slimshadyBPput$rec_percent


slimshadyBPput$odds_ratio <- unlist(lapply(1:length(slimshadyBPput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBPput$dom_count[x],(slimshadyBPput$dom_count[x]/(slimshadyBPput$dom_percent[x]/100))*(1-slimshadyBPput$dom_percent[x]/100),
                             slimshadyBPput$rec_count[x],(slimshadyBPput$rec_count[x]/(slimshadyBPput$rec_percent[x]/100))*(1-slimshadyBPput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyBPput$confint_lower <- unlist(lapply(1:length(slimshadyBPput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBPput$dom_count[x],(slimshadyBPput$dom_count[x]/(slimshadyBPput$dom_percent[x]/100))*(1-slimshadyBPput$dom_percent[x]/100),
                             slimshadyBPput$rec_count[x],(slimshadyBPput$rec_count[x]/(slimshadyBPput$rec_percent[x]/100))*(1-slimshadyBPput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyBPput$confint_higher <- unlist(lapply(1:length(slimshadyBPput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBPput$dom_count[x],(slimshadyBPput$dom_count[x]/(slimshadyBPput$dom_percent[x]/100))*(1-slimshadyBPput$dom_percent[x]/100),
                             slimshadyBPput$rec_count[x],(slimshadyBPput$rec_count[x]/(slimshadyBPput$rec_percent[x]/100))*(1-slimshadyBPput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyBPput$pval <- unlist(lapply(1:length(slimshadyBPput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyBPput$dom_count[x],(slimshadyBPput$dom_count[x]/(slimshadyBPput$dom_percent[x]/100))*(1-slimshadyBPput$dom_percent[x]/100),
                             slimshadyBPput$rec_count[x],(slimshadyBPput$rec_count[x]/(slimshadyBPput$rec_percent[x]/100))*(1-slimshadyBPput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))




slimshadyBPput$gene <- rep("gene",length(slimshadyBPput$go_term))


# taking 10 most significantly different terms
slimshadyBPput <- slimshadyBPput[order(slimshadyBPput$pval),]
slimshadyBPput <- slimshadyBPput[1:8,]
slimshadyBPput <- slimshadyBPput[order(slimshadyBPput$odds_ratio),]

#fixing abbreviated GO names
slimshadyBPput$go_term<-as.character(slimshadyBPput$go_term)
slimshadyBPput$go_term[which(slimshadyBPput$go_term=="cellular protein modification proce...")]<-"cellular protein modification process"
slimshadyBPput$go_term[which(slimshadyBPput$go_term=="cellular amino acid metabolic proce...")]<-"cellular amino acid metabolic process"
slimshadyBPput$go_term <- factor(slimshadyBPput$go_term, levels = slimshadyBPput$go_term)
#plotting heat map
d<-ggplot(data = slimshadyBPput, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45",trans="log",
                       breaks=c(min(slimshadyBPput$odds_ratio),1,max(slimshadyBPput$odds_ratio)),
                       labels=c("Enriched in PUTATIVE RECESSIVE human lethal genes","","enriched in PUTATIVE DOMINANT human lethal genes"))+  
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 20, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))






############PUTATIVE lethal genes dominant vs recessive CC################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y")&universe_df$constrained=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "CC")

slimshadyCCput <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y")&universe_df$constrained=="N"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "CC")

slimshadyCCput$rec_count <- a$Count
slimshadyCCput$rec_percent <- a$Percent
rm(a)

slimshadyCCput<-slimshadyCCput[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyCCput$dom_percent[which(slimshadyCCput$dom_count==0)]<-slimshadyCCput$dom_percent[which(slimshadyCCput$dom_count==2)][1]/2
slimshadyCCput$rec_percent[which(slimshadyCCput$rec_percent==0)]<-slimshadyCCput$rec_percent[which(slimshadyCCput$rec_count==2)][1]/2
slimshadyCCput$dom_count[which(slimshadyCCput$dom_count==0)]<-1
slimshadyCCput$rec_count[which(slimshadyCCput$rec_count==0)]<-1


slimshadyCCput$domrec_diff <- slimshadyCCput$dom_percent-slimshadyCCput$rec_percent


slimshadyCCput$odds_ratio <- unlist(lapply(1:length(slimshadyCCput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCCput$dom_count[x],(slimshadyCCput$dom_count[x]/(slimshadyCCput$dom_percent[x]/100))*(1-slimshadyCCput$dom_percent[x]/100),
                             slimshadyCCput$rec_count[x],(slimshadyCCput$rec_count[x]/(slimshadyCCput$rec_percent[x]/100))*(1-slimshadyCCput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyCCput$confint_lower <- unlist(lapply(1:length(slimshadyCCput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCCput$dom_count[x],(slimshadyCCput$dom_count[x]/(slimshadyCCput$dom_percent[x]/100))*(1-slimshadyCCput$dom_percent[x]/100),
                             slimshadyCCput$rec_count[x],(slimshadyCCput$rec_count[x]/(slimshadyCCput$rec_percent[x]/100))*(1-slimshadyCCput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyCCput$confint_higher <- unlist(lapply(1:length(slimshadyCCput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCCput$dom_count[x],(slimshadyCCput$dom_count[x]/(slimshadyCCput$dom_percent[x]/100))*(1-slimshadyCCput$dom_percent[x]/100),
                             slimshadyCCput$rec_count[x],(slimshadyCCput$rec_count[x]/(slimshadyCCput$rec_percent[x]/100))*(1-slimshadyCCput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyCCput$pval <- unlist(lapply(1:length(slimshadyCCput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyCCput$dom_count[x],(slimshadyCCput$dom_count[x]/(slimshadyCCput$dom_percent[x]/100))*(1-slimshadyCCput$dom_percent[x]/100),
                             slimshadyCCput$rec_count[x],(slimshadyCCput$rec_count[x]/(slimshadyCCput$rec_percent[x]/100))*(1-slimshadyCCput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyCCput <- slimshadyCCput[-which(slimshadyCCput$go_term=="cellular_component"),]

# taking 10 most significantly different terms
slimshadyCCput <- slimshadyCCput[order(slimshadyCCput$pval),]
slimshadyCCput <- slimshadyCCput[1:6,]
slimshadyCCput <- slimshadyCCput[order(slimshadyCCput$odds_ratio),]


slimshadyCCput$gene <- rep("gene",length(slimshadyCCput$go_term))
slimshadyCCput$go_term <- factor(slimshadyCCput$go_term, levels = slimshadyCCput$go_term)

#plotting heat map
e<-ggplot(data = slimshadyCCput, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45",trans="log",
                       breaks=c(min(slimshadyCCput$odds_ratio),1,max(slimshadyCCput$odds_ratio)),
                       labels=c("Enriched in RECESSIVE human lethal genes","","enriched in DOMINANT human lethal genes"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12),legend.position = 'none')+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))




############PUTATIVE lethal genes dominant vs recessive MF################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y")&universe_df$constrained=="Y"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "MF")

slimshadyMFput <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)

#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y")&universe_df$constrained=="N"&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "MF")

slimshadyMFput$rec_count <- a$Count
slimshadyMFput$rec_percent <- a$Percent
rm(a)

slimshadyMFput<-slimshadyMFput[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyMFput$dom_percent[which(slimshadyMFput$dom_count==0)]<-slimshadyMFput$dom_percent[which(slimshadyMFput$dom_count==1)][1]
slimshadyMFput$rec_percent[which(slimshadyMFput$rec_percent==0)]<-slimshadyMFput$rec_percent[which(slimshadyMFput$rec_count==5)][1]/5
slimshadyMFput$dom_count[which(slimshadyMFput$dom_count==0)]<-1
slimshadyMFput$rec_count[which(slimshadyMFput$rec_count==0)]<-1



slimshadyMFput$odds_ratio <- unlist(lapply(1:length(slimshadyMFput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMFput$dom_count[x],(slimshadyMFput$dom_count[x]/(slimshadyMFput$dom_percent[x]/100))*(1-slimshadyMFput$dom_percent[x]/100),
                             slimshadyMFput$rec_count[x],(slimshadyMFput$rec_count[x]/(slimshadyMFput$rec_percent[x]/100))*(1-slimshadyMFput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$estimate)
}))
slimshadyMFput$confint_lower <- unlist(lapply(1:length(slimshadyMFput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMFput$dom_count[x],(slimshadyMFput$dom_count[x]/(slimshadyMFput$dom_percent[x]/100))*(1-slimshadyMFput$dom_percent[x]/100),
                             slimshadyMFput$rec_count[x],(slimshadyMFput$rec_count[x]/(slimshadyMFput$rec_percent[x]/100))*(1-slimshadyMFput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))
slimshadyMFput$confint_higher <- unlist(lapply(1:length(slimshadyMFput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMFput$dom_count[x],(slimshadyMFput$dom_count[x]/(slimshadyMFput$dom_percent[x]/100))*(1-slimshadyMFput$dom_percent[x]/100),
                             slimshadyMFput$rec_count[x],(slimshadyMFput$rec_count[x]/(slimshadyMFput$rec_percent[x]/100))*(1-slimshadyMFput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))
slimshadyMFput$pval <- unlist(lapply(1:length(slimshadyMFput$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyMFput$dom_count[x],(slimshadyMFput$dom_count[x]/(slimshadyMFput$dom_percent[x]/100))*(1-slimshadyMFput$dom_percent[x]/100),
                             slimshadyMFput$rec_count[x],(slimshadyMFput$rec_count[x]/(slimshadyMFput$rec_percent[x]/100))*(1-slimshadyMFput$rec_percent[x]/100)),nrow=2,ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


slimshadyMFput$gene <- rep("gene",length(slimshadyMFput$go_term))
slimshadyMFput <- slimshadyMFput[-which(slimshadyMFput$go_term=="molecular_function"),]

#fixing abbreviated GO names
slimshadyMFput$go_term<-as.character(slimshadyMFput$go_term)
slimshadyMFput$go_term[which(slimshadyMFput$go_id=="GO:0016757")]<-"glycosyltransferase activity"
slimshadyMFput$go_term <- factor(slimshadyMFput$go_term, levels = slimshadyMFput$go_term)



# taking 10 most significantly different terms
slimshadyMFput <- slimshadyMFput[order(slimshadyMFput$pval),]
slimshadyMFput <- slimshadyMFput[1:8,]
slimshadyMFput <- slimshadyMFput[order(slimshadyMFput$odds_ratio),]
slimshadyMFput$go_term <- factor(slimshadyMFput$go_term, levels = slimshadyMFput$go_term)

#plotting heat map
f<-ggplot(data = slimshadyMFput, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio,width=1)) +
  scale_fill_gradient2(low = "#a020f0",mid="white",high = "#008b45",trans="log",
                       breaks=c(min(slimshadyMFput$odds_ratio),1,max(slimshadyMFput$odds_ratio)),
                       labels=c("UNCONSTRAINED","","CONSTRAINED"))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=12))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 3, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))





#saving plot####
ggarrange(ap,b,c,d,e,f,widths=c(1,1.1,1.65))
ggsave("output/Figures/5.pdf",height=30, width=26, units='cm')


rm(ap,b,c,d,e,f,slim,myCollection,slimshadyBP,slimshadyBPput,slimshadyCC,slimshadyCCput,slimshadyMF,slimshadyMFput,fl,myIds)