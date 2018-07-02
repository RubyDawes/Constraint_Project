cell_ess_hpo<- universe_df$hpo_terms[which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")]
cell_ess_genes <- universe_df$gene[which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")]
cell_uness_hpo <- universe_df$hpo_terms[which(universe_df$omim=="Y"&universe_df$cell_essential=="N")]
cell_uness_genes <- universe_df$gene[which(universe_df$omim=="Y"&universe_df$cell_essential=="N")]

hp <- get.ontology("Gene_lists/HPO/hp.obo",qualifier="HP")
library("hpoPlot")

general_cats <- c("HP:0000707","HP:0000478","HP:0000152","HP:0000119","HP:0000924","HP:0001939","HP:0003011",
                  "HP:0001871","HP:0001626","HP:0002664","HP:0002715","HP:0001574","HP:0040064","HP:0025031",
                  "HP:0000598","HP:0000818","HP:0002086","HP:0001197","HP:0003549","HP:0000769","HP:0001507","HP:0045027")
general_cats_names <- get.shortened.names(hp,general_cats)

cell_ess_hpo_anc <- lapply(cell_ess_hpo,function(x) get.ancestors(hp,x))
general_cats_celless_props <- c()
for (i in 1:length(general_cats)) {
  a <- unlist(lapply(cell_ess_hpo_anc,function(x) ifelse(general_cats[i]%in%x,1,0)))
  general_cats_celless_props[i] <- length(which(a==1))/length(a)
}

cell_uness_hpo_anc <- lapply(cell_uness_hpo,function(x) get.ancestors(hp,x))
general_cats_celluness_props <- c()
for (i in 1:length(general_cats)) {
  a <- unlist(lapply(cell_uness_hpo_anc,function(x) ifelse(general_cats[i]%in%x,1,0)))
  general_cats_celluness_props[i] <- length(which(a==1))/length(a)
}

data.frame(general_cats,general_cats_names,general_cats_celless_props,general_cats_celluness_props)
