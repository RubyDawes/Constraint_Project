lof_intol_hpo<- universe_df$hpo_terms[which(universe_df$omim=="Y"&universe_df$pLI>=0.9)]
lof_intol_genes <- universe_df$gene[which(universe_df$omim=="Y"&universe_df$pLI>=0.9)]
lof_tol_hpo <- universe_df$hpo_terms[which(universe_df$omim=="Y"&universe_df$pLI<=0.1)]
lof_tol_genes <- universe_df$gene[which(universe_df$omim=="Y"&universe_df$pLI<=0.1)]

hp <- get.ontology("Gene_lists/HPO/hp.obo",qualifier="HP")
library("hpoPlot")

general_cats <- c("HP:0000707","HP:0000478","HP:0000152","HP:0000119","HP:0000924","HP:0001939","HP:0003011",
                  "HP:0001871","HP:0001626","HP:0002664","HP:0002715","HP:0001574","HP:0040064","HP:0025031",
                  "HP:0000598","HP:0000818","HP:0002086","HP:0001197","HP:0003549","HP:0000769","HP:0001507","HP:0045027")
general_cats_names <- get.shortened.names(hp,general_cats)

lof_intol_hpo_anc <- lapply(lof_intol_hpo,function(x) get.ancestors(hp,x))
general_cats_lofintol_props <- c()
for (i in 1:length(general_cats)) {
  a <- unlist(lapply(lof_intol_hpo_anc,function(x) ifelse(general_cats[i]%in%x,1,0)))
  general_cats_lofintol_props[i] <- length(which(a==1))/length(a)
}

lof_tol_hpo_anc <- lapply(lof_tol_hpo,function(x) get.ancestors(hp,x))
general_cats_loftol_props <- c()
for (i in 1:length(general_cats)) {
  a <- unlist(lapply(lof_tol_hpo_anc,function(x) ifelse(general_cats[i]%in%x,1,0)))
  general_cats_loftol_props[i] <- length(which(a==1))/length(a)
}


data.frame(general_cats,general_cats_names,general_cats_lofintol_props,general_cats_loftol_props)
