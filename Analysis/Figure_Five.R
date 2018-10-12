############omim genes dominant vs recessive BP################
#all genes go slim annotations- human lethal DOMINANT
all <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(all)
aa<-goSlim(myCollection, slim, "BP")

slimshadyOMIMBP <- data.frame(go_id=rownames(aa),go_term=aa$Term,dom_count = aa$Count,dom_percent=aa$Percent)
rm(aa,all)
#all genes go slim annotations- human lethal RECESSIVE
myIds <- unlist(universe_df$go_terms[
  which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="XLr")&lengths(universe_df$go_terms)>0)])
myCollection <- GOCollection(myIds)
a<-goSlim(myCollection, slim, "BP")

slimshadyOMIMBP$rec_count <- a$Count
slimshadyOMIMBP$rec_percent <- a$Percent
rm(a)

slimshadyOMIMBP<-slimshadyOMIMBP[-c(38,47,69),]
#making 0 counts 1 to avoid infinity OR
slimshadyOMIMBP$dom_percent[which(slimshadyOMIMBP$dom_count==0)]<-slimshadyOMIMBP$dom_percent[which(slimshadyOMIMBP$dom_count==1)][1]
slimshadyOMIMBP$dom_count[which(slimshadyOMIMBP$dom_count==0)]<-1

slimshadyOMIMBP$domrec_diff <- slimshadyOMIMBP$dom_percent-slimshadyOMIMBP$rec_percent


slimshadyOMIMBP$odds_ratio <- unlist(lapply(1:length(slimshadyOMIMBP$go_term),function(x) {
  fish<-fisher.test(matrix(c(slimshadyOMIMBP$dom_count[x],(slimshadyOMIMBP$dom_count[x]/(slimshadyOMIMBP$dom_percent[x]/100))*(1-slimshadyOMIMBP$dom_percent[x]/100),
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

slimshadyOMIMCC$go_term<-as.character(slimshadyOMIMCC$go_term)
slimshadyOMIMCC$go_term[which(slimshadyOMIMCC$go_term=="protein-containing complex")]<- "protein containing complex"


slimshadyOMIMCC$go_term <- factor(slimshadyOMIMCC$go_term, levels = slimshadyOMIMCC$go_term)

# taking 10 most significantly different terms
slimshadyOMIMCC <- slimshadyOMIMCC[order(slimshadyOMIMCC$pval),]
slimshadyOMIMCC <- slimshadyOMIMCC[1:6,]
slimshadyOMIMCC <- slimshadyOMIMCC[order(slimshadyOMIMCC$odds_ratio),]








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

slimshadyCCput$go_term<-as.character(slimshadyCCput$go_term)
slimshadyCCput$go_term[which(slimshadyCCput$go_term=="protein-containing complex")]<- "protein containing complex"

slimshadyCCput$go_term <- factor(slimshadyCCput$go_term, levels = slimshadyCCput$go_term)





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



#plotting ####

#finding scale for legend
min<-min(c(slimshadyBPput$odds_ratio,slimshadyCCput$odds_ratio,slimshadyMFput$odds_ratio,slimshadyBP$odds_ratio,slimshadyCC$odds_ratio,slimshadyMF$odds_ratio))
max<-max(c(slimshadyBPput$odds_ratio,slimshadyCCput$odds_ratio,slimshadyMFput$odds_ratio,slimshadyBP$odds_ratio,slimshadyCC$odds_ratio,slimshadyMF$odds_ratio))

slimshadyOMIMBP$odds_ratio_adjusted<- ifelse(slimshadyOMIMBP$odds_ratio<0.1,0.1,ifelse(slimshadyOMIMBP$odds_ratio>3,3,slimshadyOMIMBP$odds_ratio))
slimshadyOMIMMF$odds_ratio_adjusted<- ifelse(slimshadyOMIMMF$odds_ratio<0.1,0.1,ifelse(slimshadyOMIMMF$odds_ratio>3,3,slimshadyOMIMMF$odds_ratio))
slimshadyOMIMCC$odds_ratio_adjusted<- ifelse(slimshadyOMIMCC$odds_ratio<0.1,0.1,ifelse(slimshadyOMIMCC$odds_ratio>3,3,slimshadyOMIMCC$odds_ratio))

slimshadyBP$odds_ratio_adjusted<- ifelse(slimshadyBP$odds_ratio<0.1,0.1,ifelse(slimshadyBP$odds_ratio>3,3,slimshadyBP$odds_ratio))
slimshadyMF$odds_ratio_adjusted<- ifelse(slimshadyMF$odds_ratio<0.1,0.1,ifelse(slimshadyMF$odds_ratio>3,3,slimshadyMF$odds_ratio))
slimshadyCC$odds_ratio_adjusted<- ifelse(slimshadyCC$odds_ratio<0.1,0.1,ifelse(slimshadyCC$odds_ratio>3,3,slimshadyCC$odds_ratio))
slimshadyBPput$odds_ratio_adjusted<- ifelse(slimshadyBPput$odds_ratio<0.1,0.1,ifelse(slimshadyBPput$odds_ratio>3,3,slimshadyBPput$odds_ratio))
slimshadyMFput$odds_ratio_adjusted<- ifelse(slimshadyMFput$odds_ratio<0.1,0.1,ifelse(slimshadyMFput$odds_ratio>3,3,slimshadyMFput$odds_ratio))
slimshadyCCput$odds_ratio_adjusted<- ifelse(slimshadyCCput$odds_ratio<0.1,0.1,ifelse(slimshadyCCput$odds_ratio>3,3,slimshadyCCput$odds_ratio))

#plotting heat map
api<-ggplot(data = slimshadyOMIMBP, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=10))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 30, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

#plotting heat map
bi<-ggplot(data = slimshadyOMIMMF, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=10),legend.position = 'none')+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))


#plotting heat map
ci<-ggplot(data = slimshadyOMIMCC, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=10))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 5, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

#plotting heat map
ap<-ggplot(data = slimshadyBP, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=10))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 30, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))


#plotting heat map
b<-ggplot(data = slimshadyMF, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=10),legend.position = 'none')+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))


#plotting heat map
c<-ggplot(data = slimshadyCC, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=10))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 5, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

  #plotting heat map
d<-ggplot(data = slimshadyBPput, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),legend.position = 'none',axis.text.y=element_text(size=10))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 30, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))


#plotting heat map
e<-ggplot(data = slimshadyMFput, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=10),legend.position = 'none')+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 25, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

#plotting heat map
f<-ggplot(data = slimshadyCCput, aes(x= gene,y = go_term)) +
  geom_tile(aes(fill = odds_ratio_adjusted,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "legend", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=10))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))





#saving plot####
ggarrange(api,bi,ci,ap,b,c,d,e,f,widths=c(1.2,1,1.1))
ggsave("output/Figures/5.pdf",height=23, width=17, units='cm')


#rm(ap,b,c,d,e,f,slim,myCollection,slimshadyBP,slimshadyBPput,slimshadyCC,slimshadyCCput,slimshadyMF,slimshadyMFput,fl,myIds)