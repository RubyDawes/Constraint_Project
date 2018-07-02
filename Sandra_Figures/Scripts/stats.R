#chi squared test on relationship between constraint and disease
M <- matrix(c(length(which(universe_df$pLI>=0.9&universe_df$omim=="Y")),
              length(which(universe_df$pLI>=0.9&is.na(universe_df$omim))),
              length(which(universe_df$pLI<0.9&universe_df$omim=="Y")),
              length(which(universe_df$pLI<0.9&is.na(universe_df$omim)))
              ),nrow=2,ncol=2)

chisq.test(M,correct=FALSE)

ks.test(universe_df$pLI[which(universe_df$omim=="Y"&!is.na(universe_df$exac))],
        universe_df$pLI[which(is.na(universe_df$omim)&!is.na(universe_df$exac))])

ggplot(universe_df,aes(factor(omim),pLI))+geom_violin(aes(fill=factor(omim)))

ggplot(universe_df,aes(factor(omim),mis_z))+geom_violin(aes(fill=factor(omim)))