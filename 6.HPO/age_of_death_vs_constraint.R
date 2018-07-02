
#how many genes causing death in infancy death are LoF intolerant

#how many genes causing Neonatal death are LoF intolerant

#how many genes causing Sudden death are LoF intolerant

#how many genes causing Death in childhood are LoF intolerant

#how many genes causing Stillbirth are LoF intolerant
length(which(grepl("Stillbirth",universe_df$hpo_names)&universe_df$pLI>=0.9))/length(which(grepl("Stillbirth",universe_df$hpo_names)&!is.na(universe_df$pLI)))

#how many genes causing Neonatal death are LoF intolerant
