
#how many AR genes show neither missense nor LoF constraint
length(which(universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))

#how many AD genes show neither missense nor LoF constraint
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$constrained=="N"))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))

#how many XL genes show neither missense nor LoF constraint
length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$constrained=="N"))/
  length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))

#how many constrained genes are OMIM genes
length(which(universe_df$mis_z>=3.09&universe_df$pLI>=0.9&universe_df$omim=="Y"))/length(which(universe_df$mis_z>=3.09&universe_df$pLI>=0.9))

#how many genes had phenotype data- 44%
length(which(universe_df$mouse_ko=="Y"))/length(universe_df$gene)
#how many oMIM genes had phenotype data-71%
length(which(universe_df$mouse_ko=="Y"&universe_df$omim=="Y"))/length(which(universe_df$omim=="Y"))
#how many non-oMIM genes had phenotype data-71%
length(which(universe_df$mouse_ko=="Y"&is.na(universe_df$omim)))/length(which(is.na(universe_df$omim)))
#how many genes were lethal in a mouse
length(which(universe_df$lethal_mouse=="Y"))/length(which(universe_df$mouse_ko=="Y"))
#how many of these were OMIM, how many were non-OMIM
length(which(universe_df$lethal_mouse=="Y"&is.na(universe_df$omim)))/length(which(is.na(universe_df$omim)&!is.na(universe_df$mouse_ko)))
length(which(universe_df$lethal_mouse=="Y"&universe_df$omim=="Y"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$mouse_ko)))

#•	ExAC: LoF – 77.4% constrained genes show no human phenotype. Of the 62.2% of these for which phenotyping information was available
#– 56.5% ‘lethal phenotype’. Implies likely to cause lethal developmental disorder in humans. Good candidates for recurrent miscarriage etc
length(which(universe_df$pLI>=0.9&is.na(universe_df$omim)))/length(which(universe_df$pLI>=0.9))
length(which(universe_df$pLI>=0.9&is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))/length(which(universe_df$pLI>=0.9&is.na(universe_df$omim)))
length(which(universe_df$pLI>=0.9&is.na(universe_df$omim)&universe_df$mouse_ko=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$pLI>=0.9&is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))

length(which(universe_df$mis_z>=3.09&is.na(universe_df$omim)))/length(which(universe_df$mis_z>=3.09))
length(which(universe_df$mis_z>=3.09&is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))/length(which(universe_df$mis_z>=3.09&is.na(universe_df$omim)))
length(which(universe_df$mis_z>=3.09&is.na(universe_df$omim)&universe_df$mouse_ko=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$mis_z>=3.09&is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))




