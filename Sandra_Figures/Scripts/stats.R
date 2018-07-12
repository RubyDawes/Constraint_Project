#chi squared test on relationship between constraint and disease
M <- matrix(c(
              length(which(universe_df$pLI>=0.9&is.na(universe_df$omim))),
              length(which(universe_df$pLI<0.9&is.na(universe_df$omim))),
              length(which(universe_df$pLI>=0.9&universe_df$omim=="Y")),
              length(which(universe_df$pLI<0.9&universe_df$omim=="Y"))
              
              ),nrow=2,ncol=2)
colnames(M) = c("Non-OMIM", "OMIM") 
rownames(M) = c("LoF constraint", "No LoF constraint") 

chisq.test(M,correct=FALSE)

M <- matrix(c(
  length(which(universe_df$constrained=="Y"&is.na(universe_df$omim))),
  length(which(universe_df$constrained=="N"&is.na(universe_df$omim))),
  length(which(universe_df$constrained=="Y"&universe_df$omim=="Y")),
  length(which(universe_df$constrained=="N"&universe_df$omim=="Y"))
  
),nrow=2,ncol=2)
colnames(M) = c("Non-OMIM", "OMIM") 
rownames(M) = c("constraint", "No constraint") 

chisq.test(M,correct=FALSE)

ks.test(universe_df$pLI[which(universe_df$omim=="Y"&!is.na(universe_df$exac))],
        universe_df$pLI[which(is.na(universe_df$omim)&!is.na(universe_df$exac))])

ggplot(universe_df,aes(factor(omim),pLI))+geom_violin(aes(fill=factor(omim)))

ggplot(universe_df,aes(factor(omim),mis_z))+geom_violin(aes(fill=factor(omim)))

N <- matrix(c(length(which(universe_df$lethal_mouse=="Y"&universe_df$omim=="Y")),
              length(which(universe_df$lethal_mouse=="Y"&is.na(universe_df$omim))),
              length(which(universe_df$lethal_mouse=="N"&universe_df$omim=="Y")),
              length(which(universe_df$lethal_mouse=="N"&is.na(universe_df$omim)))
),nrow=2,ncol=2)

chisq.test(N,correct=FALSE)

M <- matrix(c(
  3000,
  3000,
  5500,
  5000
),nrow=2,ncol=2)
colnames(M) = c("Non-OMIM", "OMIM") 
rownames(M) = c("LoF constraint", "No LoF constraint") 

chisq.test(M,correct=FALSE)

library(binom)
binom.confint(x=length(which(universe_df$pLI>=0.9&universe_df$omim=="Y")),n=length(which(universe_df$omim=="Y")),method='wilson')$lower

binom.confint(x=length(which(universe_df$lethal_mouse=="Y"&universe_df$omim=="Y")),n=length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse))))

binom.confint(x=2275,n=2275)

binom.confint(x=19,n=20,conf.level=0.95)

